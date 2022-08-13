# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 18:31:08 2022

@author: rikci
"""

import numpy as np
from math import pi, sqrt, acos
from qiskit import QuantumCircuit, transpile, IBMQ
from qiskit.providers.aer import AerSimulator, QasmSimulator
import qiskit.providers.aer.noise as noise
from qiskit.circuit.library import QFT
from qiskit.visualization import plot_histogram
from qiskit.tools import job_monitor
from IPython.display import display
from lattice import Lattice
import phases

#IBMQ.load_account()
provider = IBMQ.get_provider(hub='ibm-q')
backend = provider.get_backend('ibm_nairobi')
#noise_model = noise.NoiseModel.from_backend(backend)

# Error probabilities
prob_1 = 0.003  # 1-qubit gate
prob_2 = 0.02   # 2-qubit gate
prob_3 = 1 - ((1-prob_2)**13)*((1-prob_1)**12) # 3-qubit gate

# Depolarizing quantum errors
error_1 = noise.depolarizing_error(prob_1, 1)
error_2 = noise.depolarizing_error(prob_2, 2)
error_3 = noise.depolarizing_error(prob_3, 3)

# Add errors to noise model
noise_model = noise.NoiseModel()
noise_model.add_all_qubit_quantum_error(error_1, ['rz', 'sx', 'x', 'h', 'z', 's', 'sdg'])
noise_model.add_all_qubit_quantum_error(error_2, ['cx', 'cp'])
noise_model.add_all_qubit_quantum_error(error_3, ['ccx'])
print(noise_model)

"""
IMPORTANT: the encoding between links of the lattice and qubits is, for Z3:
    QUBIT 0 = 2*LATTICE
    QUBIT 1 = 2*LATTICE + 1
Except for the ancilla, which is the qubit N° N_qubits-1
"""


class ToricCode(Lattice):
    
    def __init__(self, size, pbc, spl, ancillas):
        Lattice.__init__(self, size, pbc, spl)
        self.N_ancillas = ancillas
        self.N_qubits = 2*self.nlinks + self.N_ancillas
        self.N_shots = 2000
        self.circuit = QuantumCircuit(self.N_qubits, ancillas)
        
    def qubits_from_links(self, link_list):
        """Returns the qubits associated to a given list of links"""
        qubit_list = [(2*link, 2*link+1) for link in link_list]
        return qubit_list

    def Z(self, qubit_pair, power = 1, control_qubit = -1, qc = None):
        if power%3 == 0:
            return
        if qc is None:
            qc = self.circuit
        if control_qubit == -1:
            if power%3 == 1:
                qc.x(qubit_pair[0])
                qc.cx(qubit_pair[0], qubit_pair[1])
                qc.cx(qubit_pair[1], qubit_pair[0])
            elif power%3 == 2:
                qc.cx(qubit_pair[1], qubit_pair[0])
                qc.cx(qubit_pair[0], qubit_pair[1])
                qc.x(qubit_pair[0])
        else:
            if power%3 == 1:
                qc.cx(control_qubit, qubit_pair[0])
                qc.ccx(control_qubit, qubit_pair[0], qubit_pair[1])
                qc.ccx(control_qubit, qubit_pair[1], qubit_pair[0])
            elif power%3 == 2:
                qc.ccx(control_qubit, qubit_pair[1], qubit_pair[0])
                qc.ccx(control_qubit, qubit_pair[0], qubit_pair[1])
                qc.cx(control_qubit, qubit_pair[0])
        
    def X(self, qubit_pair, power = 1, control_qubit = -1, qc = None):
        if power%3 == 0:
            return
        if qc is None:
            qc = self.circuit
        if control_qubit == -1:
            if power%3 == 1:
                qc.p(-2*pi/3, qubit_pair[0])
                qc.p(2*pi/3, qubit_pair[1])
            elif power%3 == 2:
                qc.p(2*pi/3, qubit_pair[0])
                qc.p(-2*pi/3, qubit_pair[1])
        else:
            if power%3 == 1:
                qc.cp(-2*pi/3, control_qubit, qubit_pair[0])
                qc.cp(2*pi/3, control_qubit, qubit_pair[1])
            elif power%3 == 2:
                qc.cp(2*pi/3, control_qubit, qubit_pair[0])
                qc.cp(-2*pi/3, control_qubit, qubit_pair[1])
            
            
    
    def MagneticGS(self, qc = None):
        if qc is None:
            qc = self.circuit
        
        #all rows in parallel, three steps, sequential in the columns
        for x in range(self.Lx - 1):
            for i in range(3):
                #three iterations to perform in parallel these steps on each column but the last
                #i=0: fourier transform the right link and act on right->bottom links
                #i=1: act on right->top links
                #i=2: act on right->left links (this is the only step that must
                #be done sequentially in the columns)
                for y in range(self.Ly):
                    link_list = self.plaquette((x,y), from_zero = True)
                    qubit_list = self.qubits_from_links(link_list)
                    fourier_pair = qubit_list[1]
                    pair_index = i+1
                    
                    #if first step:
                    if i==0:    
                        qc.ry(2*acos(1/sqrt(3)), fourier_pair[0])
                        qc.ch(fourier_pair[0], fourier_pair[1])
                        qc.x(fourier_pair[0])
                        pair_index = 0
                    
                    pair = qubit_list[pair_index]
                    qc.cx(fourier_pair[0], fourier_pair[1])
                    if pair_index < 2:
                        self.Z(pair, power = 1, control_qubit = fourier_pair[0])
                        self.Z(pair, power = 1, control_qubit = fourier_pair[1])
                    else:
                        self.Z(pair, power = 2, control_qubit = fourier_pair[1])
                        self.Z(pair, power = 2, control_qubit = fourier_pair[0])
                    qc.cx(fourier_pair[0], fourier_pair[1])
                                          
        #last column from top to bottom, last plaquette excluded       
        x = self.Lx - 1
        for y in range(self.Ly - 1, 0, -1):
            link_list = self.plaquette((x,y), from_zero = True)
            qubit_list = self.qubits_from_links(link_list)
            fourier_pair = qubit_list[0]
            
            #fourier transform the lower link
            qc.ry(2*acos(1/sqrt(3)), fourier_pair[0])
            qc.ch(fourier_pair[0], fourier_pair[1])
            qc.x(fourier_pair[0])

            for pair_index, pair in enumerate(qubit_list):
                if pair_index > 0:
                    qc.cx(fourier_pair[0], fourier_pair[1])
                    if pair_index < 2:
                        self.Z(pair, power = 1, control_qubit = fourier_pair[0])
                        self.Z(pair, power = 1, control_qubit = fourier_pair[1])
                    else:
                        self.Z(pair, power = 2, control_qubit = fourier_pair[1])
                        self.Z(pair, power = 2, control_qubit = fourier_pair[0])
                    qc.cx(fourier_pair[0], fourier_pair[1])
                    
                    
    def Z_string(self, site0, site1, power = 1, control_qubit = -1, circuit = None):        
        if power%3==0:
            return
        delta = site1[0] - site0[0] + site1[1] - site0[1]
        if delta==0:
            raise RuntimeError("Z string of lenght zero, please insert valid path")
        elif delta > 0:
            path = self.path(site0, site1)
        else:
            path = self.path(site1, site0)
            power = 3 - power%3
        qubit_path = self.qubits_from_links(path)
        for qubit_pair in qubit_path:
            self.Z(qubit_pair, power, control_qubit, circuit)    
                    
    def X_string(self, site0, site1, power = 1, control_qubit = -1, circuit = None):
        if power%3==0:
            return
        delta = site1[0] - site0[0] + site1[1] - site0[1]
        if delta==0:
            raise RuntimeError("X string of lenght zero, please insert valid path")
        elif delta > 0:
            path = self.path(site0, site1)
        else:
            path = self.path(site1, site0)
            power = 3 - power%3
        qubit_path = self.qubits_from_links(path)
        for qubit_pair in qubit_path:
            self.X(qubit_pair, power, control_qubit, circuit)
    
    
    def Bp(self, site, power = 1, control_qubit = -1, circuit = None):   
        if type(power) is not int:
            raise RuntimeError("Plaquette raised to a non-integer power")
        x, y = site
        sites = [(x, y), (x+1, y), (x+1, y+1), (x, y+1)]
        if power%3 == 0:
            return
        else:
            for i in range(4):
                self.Z_string(sites[i], sites[(i+1)%4], power, control_qubit, circuit)
        
    def init_dyons(self, low_charges = (1,1), up_charges = (1,1), **kwargs):
        """      
        Parameters
        ----------
        low_charges : tuple(int, int), optional
            Charges of the lower dyon, created on the (1,1) plaquette (electric
            charge on the top-left vertex).
            The first integer is the electric charge, the second is the 
            magnetic charge. The default is (1,1).
        up_charges : tuple(int, int), optional
            Charges of the upper dyon, created on the (0,0) plaquette.
            The default is (1,1).

        Returns
        -------
        None.
        """     
        if self.Lx < 3 or self.Ly < 2:
            raise RuntimeError("Lattice too short to initialize dyons. Make sure the lattice is 3x2 at least.")
        self.low_dyon_pow = low_charges
        self.up_dyon_pow = up_charges
        #Initialize 1st, lower dyon
        self.Z_string((0,1), (2,1), power = low_charges[0], **kwargs)
        self.X_string((1,0), (1,1), power = low_charges[1], **kwargs)
        self.X_string((2,0), (2,1), power = low_charges[1], **kwargs)
        #Initialize 2nd, higher dyon
        self.Z_string((1,2), (2,2), power = up_charges[0], **kwargs)
        self.X_string((2,1), (2,2), power = up_charges[1], **kwargs)
        
        
    def exchange_countclock(self, N_exchanges = 1, **kwargs):
        if self.low_dyon_pow != self.up_dyon_pow:
            print("\nWarning: the dyons are not identical. Only full rotations are meaningful.")
        for i in range(N_exchanges):
            low = self.low_dyon_pow
            up = self.up_dyon_pow
            #1st to the right
            self.Z_string((1,1), (0,1), power = low[0], **kwargs)
            self.X_string((1,1), (1,0), power = low[1], **kwargs)
            #2nd to the left and down
            self.Z_string((0,2), (1,2), power = up[0], **kwargs)
            self.Z_string((0,1), (0,2), power = up[0], **kwargs)
            self.X_string((1,1), (1,2), power = up[1], **kwargs)
            self.X_string((1,1), (0,1), power = up[1], **kwargs)
            #1st up, complete exchange
            self.Z_string((1,2), (1,1), power = low[0], **kwargs)
            self.X_string((1,1), (2,1), power = low[1], **kwargs)
            #update positions
            self.low_dyon_pow = up
            self.up_dyon_pow = low
    

test = ToricCode((3,2), (True,True), 2, 3)
test.MagneticGS()
ancilla = test.N_qubits-test.N_ancillas
reg_list = [ancilla, ancilla+1, ancilla+2]
test.circuit.h(reg_list)

#Initialize e and m particles and rotate them

test.Z_string((1,1), (2,1))
test.X_string((1,0), (1,1), power = 1)
test.X_string((2,0), (2,1), power = 1)
test.Bp((0,0), power=2, control_qubit=reg_list[0])
test.Bp((0,0), power=4, control_qubit=reg_list[1])
test.Bp((0,0), power=8, control_qubit=reg_list[2])


#Access non trivial sector and measure it with 't Hooft loop
"""
test.Z_string((0,1), (3,1), power = 2)
test.X_string((1,0), (2,0), control_qubit=ancilla)
test.X_string((1,1), (2,1), control_qubit=ancilla)
"""
#Initialize dyons and exchange them an arbitrary number of times
"""
test.init_dyons(low_charges = (1,1), up_charges = (1,1))
test.exchange_countclock(2, control_qubit = ancilla)
"""


#test.circuit.h(reg_list)
test.circuit.append(QFT(3, inverse=True), [test.circuit.qubits[reg] for reg in reg_list])
test.circuit.measure(reg_list, [0,1,2])
qc_list = [test.circuit]

"""
circuit2 = test.circuit.copy()
circuit2.sdg(ancilla)
qc_list = [test.circuit, circuit2]
for circuit in qc_list:
    circuit.h(ancilla)
    circuit.measure(ancilla, ancilla)
"""

#Draw the circuit and run it through the simulator with a given number of shots
display(test.circuit.draw('mpl'))
#simulator = QasmSimulator()
simulator = AerSimulator(method = "statevector")
tqc_list = transpile(qc_list, simulator)

print("\nTest circuit depth:", qc_list[0].depth())
print("Test circuit size:", qc_list[0].size())
print("Test circuit width:", qc_list[0].width())
print("\nTranspiled circuit depth:", tqc_list[0].depth())
print("Transpiled circuit size:", tqc_list[0].size())
print("Transpiled circuit width:", tqc_list[0].width())


job = simulator.run(tqc_list, job_name = "Toric Test", shots = test.N_shots, memory = True, noise_model = noise_model)
job_monitor(job)
result = job.result()

counts = result.get_counts()
print("\nTotal counts are:", counts)
display(plot_histogram(counts))

print("\nRunning time {}s".format(result.time_taken))


    

