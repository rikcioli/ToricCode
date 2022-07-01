# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 18:31:08 2022

@author: rikci
"""

import numpy as np
from math import pi
from qiskit import QuantumCircuit, transpile
from qiskit.providers.aer import AerSimulator
from IPython.display import display
from lattice import Lattice

"""
IMPORTANT: the encoding between links of the lattice and qubits is, for Z4:
    QUBIT 0 = 2*LATTICE
    QUBIT 1 = 2*LATTICE + 1
Except for the ancilla, which is the qubit N° N_qubits-1
"""


class ToricCode(Lattice):
    def __init__(self, size, pbc, spl, ancillas):
        Lattice.__init__(self, size, pbc, spl)
        self.N_ancillas = ancillas
        self.N_qubits = 2*self.nlinks + self.N_ancillas
        self.N_shots = 10000
        self.circuit = QuantumCircuit(self.N_qubits, self.N_qubits)
        
    def qubits_from_links(self, link_list):
        """Returns the qubits associated to a given list of links"""
        qubit_list = [(2*link, 2*link+1) for link in link_list]
        return qubit_list
        
    def ElectricGS(self):
        for link in range(self.nlinks):
            qubit_list = self.qubits_from_links([link])
            for pair in qubit_list:
                self.circuit.h(pair[0])
                self.circuit.h(pair[1])               
    
    def MagneticGS(self):
        """
        Initializes magnetic ground state, by applying (1+Bp) on each plaquette.
        The application proceeds from bottom to top, using the upper link as target of
        the other 3, until the upper row is reached.
        It then proceeds from left to right, using the right link as target.
        The last plaquette is untouched, since already determined by the PBC.
        All the non-target qubits need to be initialized to eigenstates
        of X with eigenvalue +1.    
        """
        self.ElectricGS()
        #all rows in parallel, three steps, sequential in the columns
        for x in range(self.Lx - 1):
            for i in range(3):
                for y in range(self.Ly):
                    link_list = self.plaquette((x,y), from_zero = True)
                    qubit_list = self.qubits_from_links(link_list)
                    pair_index = i+1
                    #if first step, hadamard the right link and act on bottom->right links
                    if i==0:    
                        self.circuit.h(qubit_list[1][0])
                        self.circuit.h(qubit_list[1][1])
                        pair_index = 0
                        
                    pair = qubit_list[pair_index]
                    self.circuit.h(pair[1])                         #H(2)
                    if pair_index > 1: self.circuit.x(pair[1])      #X(2) if up or left link
                    self.circuit.ccx(pair[1], pair[0], qubit_list[1][1]) #CCX(2,1,4)
                    if pair_index > 1: self.circuit.x(pair[1])      #X(2) if up or left link
                    self.circuit.h(pair[1])                         #H(2)
                    self.circuit.cx(pair[0], qubit_list[1][0])      #CX(1,3)
                    self.circuit.cx(pair[1], qubit_list[1][1])      #CX(2,4)
                    
                        
        #last column from top to bottom, last plaquette excluded
        
        x = self.Lx - 1
        for y in range(self.Ly - 1, 0, -1):
            link_list = self.plaquette((x,y), from_zero = True)
            qubit_list = self.qubits_from_links(link_list)
            self.circuit.h(qubit_list[0][0])
            self.circuit.h(qubit_list[0][1])
            for pair_index, pair in enumerate(qubit_list):
                if pair_index > 0:
                    self.circuit.h(pair[1])                  #H(2)
                    if pair_index > 1: self.circuit.x(pair[1])      #X(2) if up or left link
                    self.circuit.ccx(pair[1], pair[0], qubit_list[0][1]) #CCX(2,1,4)
                    if pair_index > 1: self.circuit.x(pair[1])      #X(2) if up or left link
                    self.circuit.h(pair[1])                  #H(2)
                    self.circuit.cx(pair[0], qubit_list[0][0])   #CX(1,3)
                    self.circuit.cx(pair[1], qubit_list[0][1])   #CX(2,4)
        

            
    def CZ_oldstring(self, site0, site1, control_qubit = 0):
        path = self.path(site0, site1)
        delta = site1[0] - site0[0] + site1[1] - site0[1]
        if delta==0:
            raise RuntimeError("Z string of lenght zero, please insert valid path")
        qubit_path = self.qubits_from_links(path)
        if control_qubit == 0:        
            if delta > 0:
                for qubit in qubit_path:
                    self.circuit.z(qubit[0])
                    self.circuit.s(qubit[1])
            else:
                for qubit in qubit_path:
                    self.circuit.z(qubit[0])
                    self.circuit.sdg(qubit[1])
        else:
            if delta > 0:
                for qubit in qubit_path:
                    self.circuit.cz(control_qubit, qubit[0])
                    self.circuit.cp(pi/2, control_qubit, qubit[1])                    
            else:
                for qubit in qubit_path:
                    self.circuit.cz(control_qubit, qubit[0])
                    self.circuit.cp(-pi/2, control_qubit, qubit[1])
                    
    def CX_string(self, site0, site1, control_qubit = 0):
        path = self.path(site0, site1)
        delta = site1[0] - site0[0] + site1[1] - site0[1]
        if delta==0:
            raise RuntimeError("Z string of lenght zero, please insert valid path")
        qubit_path = self.qubits_from_links(path)
        if control_qubit == 0:        
            if delta > 0:
                for qubit in qubit_path:
                    self.circuit.x(qubit[0])
                    self.circuit.h(qubit[1])
                    self.circuit.s(qubit[1])
                    self.circuit.h(qubit[1])
            else:
                for qubit in qubit_path:
                    self.circuit.x(qubit[0])
                    self.circuit.h(qubit[1])
                    self.circuit.sdg(qubit[1])
                    self.circuit.h(qubit[1])
        else:
            if delta > 0:
                for qubit in qubit_path:                  
                    self.circuit.cx(control_qubit, qubit[0])
                    self.circuit.h(qubit[1])
                    self.circuit.cp(pi/2, control_qubit, qubit[1]) 
                    self.circuit.h(qubit[1])
            else:
                for qubit in qubit_path:
                    self.circuit.cx(control_qubit, qubit[0])
                    self.circuit.h(qubit[1])
                    self.circuit.cp(-pi/2, control_qubit, qubit[1])
                    self.circuit.h(qubit[1])
            
    def CX_oldstring(self, site0, site1, control_qubit = 0):
        path = self.path(site0, site1)
        delta = site1[0] - site0[0] + site1[1] - site0[1]
        if delta==0:
            raise RuntimeError("Z string of lenght zero, please insert valid path")
        qubit_path = self.qubits_from_links(path)
        if control_qubit == 0:        
            if delta > 0:
                for qubit in qubit_path:
                    self.circuit.cx(qubit[1], qubit[0])
                    self.circuit.x(qubit[0])
                    self.circuit.x(qubit[1])
            else:
                for qubit in qubit_path:
                    self.circuit.x(qubit[0])
                    self.circuit.x(qubit[1])
                    self.circuit.cx(qubit[1], qubit[0])
        else:
            if delta > 0:
                for qubit in qubit_path:
                    self.circuit.ccx(control_qubit, qubit[1], qubit[0])
                    self.circuit.cx(control_qubit, qubit[0])
                    self.circuit.cx(control_qubit, qubit[1])
                    
            else:
                for qubit in qubit_path:
                    self.circuit.cx(control_qubit, qubit[0])
                    self.circuit.cx(control_qubit, qubit[1])
                    self.circuit.ccx(control_qubit, qubit[1], qubit[0])
                    
    def CZ_string(self, site0, site1, control_qubit = 0):
        path = self.path(site0, site1)
        delta = site1[0] - site0[0] + site1[1] - site0[1]
        if delta==0:
            raise RuntimeError("Z string of lenght zero, please insert valid path")
        qubit_path = self.qubits_from_links(path)
        if control_qubit == 0:        
            if delta > 0:
                for qubit in qubit_path:
                    self.circuit.z(qubit[1])
                    self.circuit.z(qubit[0])
                    self.circuit.cx(qubit[0], qubit[1])
                    
            else:
                for qubit in qubit_path:
                    self.circuit.cx(qubit[0], qubit[1])
                    self.circuit.z(qubit[0])
                    self.circuit.z(qubit[1])
        else:
            if delta > 0:
                for qubit in qubit_path:
                    self.circuit.cz(control_qubit, qubit[1])
                    self.circuit.cz(control_qubit, qubit[0])
                    self.circuit.ccx(control_qubit, qubit[0], qubit[1])
                                       
            else:
                for qubit in qubit_path:
                    self.circuit.ccx(control_qubit, qubit[0], qubit[1])
                    self.circuit.cz(control_qubit, qubit[0])
                    self.circuit.cz(control_qubit, qubit[1])
                    
                    
    def CZsq_string(self, site0, site1, control_qubit = 0):
        path = self.path(site0, site1)
        delta = site1[0] - site0[0] + site1[1] - site0[1]
        if delta==0:
            raise RuntimeError("Z string of lenght zero, please insert valid path")
        qubit_path = self.qubits_from_links(path)
        if control_qubit == 0:                    
            for qubit in qubit_path:
                self.circuit.z(qubit[0])                    
        else:
            for qubit in qubit_path:
                self.circuit.cz(control_qubit, qubit[0])
    
    def CXsq_string(self, site0, site1, control_qubit = 0):
        path = self.path(site0, site1)
        delta = site1[0] - site0[0] + site1[1] - site0[1]
        if delta==0:
            raise RuntimeError("Z string of lenght zero, please insert valid path")
        qubit_path = self.qubits_from_links(path)
        if control_qubit == 0:                    
            for qubit in qubit_path:
                self.circuit.x(qubit[1])                    
        else:
            for qubit in qubit_path:
                self.circuit.cx(control_qubit, qubit[1])
    
    def Bp(self, site, power = 1, control_qubit = 0):   
        if type(power) is not int:
            raise RuntimeError("Plaquette raised to a non-integer power")
        x, y = site
        sites = [(x, y), (x+1, y), (x+1, y+1), (x, y+1)]
        if power%4 == 0:
            return
        elif power%4 == 1:            
            for i in range(4):
                self.CZ_string(sites[i], sites[(i+1)%4], control_qubit)
        elif power%4 == 2:
            for i in range(4):
                self.CZsq_string(sites[i], sites[(i+1)%4], control_qubit)
        else:
            sites = [sites[0], sites[3], sites[2], sites[1]]
            for i in range(4):
                self.CZ_string(sites[i], sites[(i+1)%4], control_qubit)                                
                    
                    
def evaluate_phase(counts):   
    #Given the counts of a measurement of the ancilla, returns the relative phase
    #between the 2 Z eigenstates   
    values = list(counts.values())
    keys = list(counts.keys())
    phase = 0
    if len(keys) == 1:
        if int(keys[0]) > 0:
            phase = np.pi
    else:
        cos_phi = (values[0] - values[1])/test.N_shots
        phase = np.arccos(cos_phi)
    return phase

def evaluate_sin(counts):
    values = list(counts.values())
    keys = list(counts.keys())
    sin_phi = 1
    if len(keys) == 1:
        if int(keys[0]) > 0:
            sin_phi = -1
    else:
        sin_phi = (values[0] - values[1])/test.N_shots
    return sin_phi
    

sine_eval = False

test = ToricCode((3,3), (True,True), 2, 1)
test.MagneticGS()
ancilla = test.N_qubits-test.N_ancillas

test.circuit.h(ancilla)

"""
#Initialize 1st, lower dyon
test.CZ_string((0,1), (2,1))
test.CX_string((1,0), (1,1))
test.CX_string((2,0), (2,1))

#Initialize 2nd, higher dyon
test.CZ_string((1,2), (1,3))
test.CX_string((2,2), (1,2))

#1st to the right
test.CZ_string((1,1), (0,1), ancilla)
test.CX_string((1,1), (1,0), ancilla)

#2nd to the left and down
test.CZ_string((0,2), (1,2), ancilla)
test.CX_string((1,1), (1,2), ancilla)
test.CZ_string((0,1), (0,2), ancilla)
test.CX_string((1,1), (0,1), ancilla)

#1st up, complete exchange
test.CZ_string((1,2), (1,1), ancilla)
test.CX_string((1,1), (2,1), ancilla)
"""


test.CZ_string((1,1), (2,1))
test.CX_string((1,0), (1,1))
test.CX_string((2,0), (2,1))
#test.CZ_string((0,2), (3,2))
test.Bp((2,0), power=3, control_qubit=ancilla)


if sine_eval == True: test.circuit.sdg(ancilla)
test.circuit.h(ancilla)

test.circuit.measure(10, 10)


#Draw the circuit and run it through the simulator with a given number of shots
display(test.circuit.draw())
simulator = AerSimulator(method = "statevector")
compiled_circuit = transpile(test.circuit, simulator)
job = simulator.run(compiled_circuit, shots = test.N_shots)
result = job.result()
counts = result.get_counts(compiled_circuit)
print("\nTotal counts are:",counts)
print("\nRunning time {}s".format(result.time_taken))

#Extract relative phase from the counts
if sine_eval == True:
    sin = evaluate_sin(counts)
    print("\nSine of phase is:", sin)
else:
    phase = evaluate_phase(counts)
    print("\nTotal phase is:", phase)
    

