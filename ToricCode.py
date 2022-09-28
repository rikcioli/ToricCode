# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 18:31:08 2022

@author: rikci
"""

from math import pi, sqrt, acos
from qiskit import QuantumCircuit
import qiskit.extensions as ext
from lattice import Lattice


"""
IMPORTANT: the encoding between links of the lattice and qubits is, for Z3 and Z4:
    QUBIT 0 = 2*LATTICE
    QUBIT 1 = 2*LATTICE + 1
Except for the ancilla, which is the qubit N° N_qubits-1
"""


class ToricCode(Lattice):
    
    def __init__(self, size, pbc, spl, N_ancillas = 0):
        Lattice.__init__(self, size, pbc, spl)
        self.N_ancillas = N_ancillas
        self.N_qubits = 2*self.nlinks + self.N_ancillas
        self.circuit = QuantumCircuit(self.N_qubits, N_ancillas)

        
    def qubits_from_links(self, link_list):
        """Returns the qubits associated to a given list of links"""
        qubit_list = [(2*link, 2*link+1) for link in link_list]
        return qubit_list
    
    def gates(self, noisy = True):
        """Initialize noiseless gates to use in simulations"""
        h_label = x_label = ry_label = ch_label = cx_label = ccx_label = None
        if noisy is False:
            x_label = 'X_nl'
            h_label = 'H_nl'
            ry_label = 'Ry_nl'
            ch_label = 'CH_nl'
            cx_label = ccx_label = 'nl'
        gates = {
            'x': ext.XGate(label = x_label),
            'h': ext.HGate(label = h_label),
            'ry': ext.RYGate(theta = 2*acos(1/sqrt(3)), label = ry_label),
            'cx': ext.CXGate(label = cx_label),
            'ch': ext.CHGate(label = ch_label),
            'ccx': ext.CCXGate(label = ccx_label),
            }
        return gates
        

    def Z(self, qubit_pair, power = 1, control_qubit = -1, qc = None):
        return
        
    def X(self, qubit_pair, power = 1, control_qubit = -1, qc = None):
        return         
            
    def MagneticGS(self, qc = None, noisy = True):
        return              
    
                    
    def Z_string(self, site0, site1, power = 1, control_qubit = -1, qc = None, noisy = True):        
        if power%self.spl==0:
            return
        delta = site1[0] - site0[0] + site1[1] - site0[1]
        if delta==0:
            raise RuntimeError("Z string of lenght zero, please insert valid path")
        elif delta > 0:
            path = self.path(site0, site1)
        else:
            path = self.path(site1, site0)
            power = self.spl - power%self.spl
        qubit_path = self.qubits_from_links(path)
        for qubit_pair in qubit_path:
            self.Z(qubit_pair, power, control_qubit, qc, noisy)    
                    
    def X_string(self, site0, site1, power = 1, control_qubit = -1, qc = None):
        if power%self.spl==0:
            return
        delta = site1[0] - site0[0] + site1[1] - site0[1]
        if delta==0:
            raise RuntimeError("X string of lenght zero, please insert valid path")
        elif delta > 0:
            path = self.path(site0, site1)
        else:
            path = self.path(site1, site0)
            power = self.spl - power%self.spl
        qubit_path = self.qubits_from_links(path)
        for qubit_pair in qubit_path:
            self.X(qubit_pair, power, control_qubit, qc)
     
    def Bp(self, site, power = 1, control_qubit = -1, qc = None, noisy = True):   
        if type(power) is not int:
            raise RuntimeError("Plaquette raised to a non-integer power")
        x, y = site
        sites = [(x, y), (x+1, y), (x+1, y+1), (x, y+1)]
        if power%self.spl == 0:
            return
        else:
            for i in range(4):
                self.Z_string(sites[i], sites[(i+1)%4], power, control_qubit, qc, noisy)
        
    def init_dyons(self, low_charges = (1,1), up_charges = (1,1), **kwargs):
        """      
        Parameters
        ----------
        low_charges : tuple(int, int), optional
            Charges of the lower dyon, created on the (0,0) plaquette (electric
            charge on the top-left vertex).
            The first integer is the electric charge, the second is the 
            magnetic charge. The default is (1,1).
        up_charges : tuple(int, int), optional
            Charges of the upper dyon, created on the (1,1) plaquette.
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
            if N_exchanges%2 == 1:
                raise RuntimeError("Dyons are not identical. Only an even number of exchanges is meaningful.")
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
            
            #New complete exchange, DO NOT DELETE
            #self.Z_string((2,1),(1,1), power = low[0], **kwargs)    #move e right
            #self.X_string((1,1), (2,1), power = low[1], **kwargs)   #move m up
            #self.Z_string((2,2), (2,1), power = low[0], **kwargs)   #move e up and left to complete
            #self.Z_string((1,2), (2,2), power = low[0], **kwargs)
            
            #Old complete exchange, DO NOT DELETE
            self.Z_string((1,2), (1,1), power = low[0], **kwargs)
            self.X_string((1,1), (2,1), power = low[1], **kwargs)
            
            #update positions
            self.low_dyon_pow = up
            self.up_dyon_pow = low
            

class Z3(ToricCode):
    
    def __init__(self, size, pbc, N_ancillas):
        ToricCode.__init__(self, size, pbc, 3, N_ancillas)
    
    def Z(self, qubit_pair, power = 1, control_qubit = -1, qc = None, noisy = True):
        if power%3 == 0:
            return
        if qc is None:
            qc = self.circuit
        gates = self.gates(noisy)
        if control_qubit == -1:
            if power%3 == 1:
                qc.append(gates['x'], [qubit_pair[0]])
                qc.append(gates['cx'], [qubit_pair[0], qubit_pair[1]])
                qc.append(gates['cx'], [qubit_pair[1], qubit_pair[0]])
            elif power%3 == 2:
                qc.append(gates['cx'], [qubit_pair[1], qubit_pair[0]])
                qc.append(gates['cx'], [qubit_pair[0], qubit_pair[1]])
                qc.append(gates['x'], [qubit_pair[0]])
        else:
            if power%3 == 1:
                qc.append(gates['cx'], [control_qubit, qubit_pair[0]])
                qc.append(gates['ccx'], [control_qubit, qubit_pair[0], qubit_pair[1]])
                qc.append(gates['ccx'], [control_qubit, qubit_pair[1], qubit_pair[0]])
            elif power%3 == 2:
                qc.append(gates['ccx'], [control_qubit, qubit_pair[1], qubit_pair[0]])
                qc.append(gates['ccx'], [control_qubit, qubit_pair[0], qubit_pair[1]])
                qc.append(gates['cx'], [control_qubit, qubit_pair[0]])
    
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
                
                
    def MagneticGS(self, qc = None, noisy = True):
        if qc is None:
            qc = self.circuit
        gates = self.gates(noisy)
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
                        qc.append(gates['ry'], [fourier_pair[0]])
                        qc.append(gates['ch'], [fourier_pair[0], fourier_pair[1]])
                        qc.append(gates['x'], [fourier_pair[0]])
                        pair_index = 0
                    
                    pair = qubit_list[pair_index]
                    qc.append(gates['cx'], [fourier_pair[0], fourier_pair[1]])
                    if pair_index < 2:
                        self.Z(pair, power = 1, control_qubit = fourier_pair[0], noisy=noisy)
                        self.Z(pair, power = 1, control_qubit = fourier_pair[1], noisy=noisy)
                    else:
                        self.Z(pair, power = 2, control_qubit = fourier_pair[1], noisy=noisy)
                        self.Z(pair, power = 2, control_qubit = fourier_pair[0], noisy=noisy)
                    qc.append(gates['cx'], [fourier_pair[0], fourier_pair[1]])
                                          
        #last column from top to bottom, last plaquette excluded       
        x = self.Lx - 1
        for y in range(self.Ly - 1, 0, -1):
            link_list = self.plaquette((x,y), from_zero = True)
            qubit_list = self.qubits_from_links(link_list)
            fourier_pair = qubit_list[0]
            
            #fourier transform the lower link
            qc.append(gates['ry'], [fourier_pair[0]])
            qc.append(gates['ch'], [fourier_pair[0], fourier_pair[1]])
            qc.append(gates['x'], [fourier_pair[0]])

            for pair_index, pair in enumerate(qubit_list):
                if pair_index > 0:
                    qc.append(gates['cx'], [fourier_pair[0], fourier_pair[1]])
                    if pair_index < 2:
                        self.Z(pair, power = 1, control_qubit = fourier_pair[0], noisy=noisy)
                        self.Z(pair, power = 1, control_qubit = fourier_pair[1], noisy=noisy)
                    else:
                        self.Z(pair, power = 2, control_qubit = fourier_pair[1], noisy=noisy)
                        self.Z(pair, power = 2, control_qubit = fourier_pair[0], noisy=noisy)
                    qc.append(gates['cx'], [fourier_pair[0], fourier_pair[1]])
        return
                
                
class Z4(ToricCode):
    
    def __init__(self, size, pbc, N_ancillas):
        ToricCode.__init__(self, size, pbc, 4, N_ancillas)
        
    def Z(self, qubit_pair, power = 1, control_qubit = -1, qc = None, noisy = True):
        if power%4 == 0:
            return
        if qc is None:
            qc = self.circuit
        gates = self.gates(noisy)
        if control_qubit == -1:        
            if power%4 == 2:                #Z^2 
                qc.append(gates['x'], [qubit_pair[0]])
            elif power%4 == 1:              #Z
                qc.append(gates['x'], [qubit_pair[1]])
                qc.append(gates['x'], [qubit_pair[0]])
                qc.append(gates['cx'], [qubit_pair[1], qubit_pair[0]])
            elif power%4 == 3:              #Zdag
                qc.append(gates['cx'], [qubit_pair[1], qubit_pair[0]])
                qc.append(gates['x'], [qubit_pair[0]])
                qc.append(gates['x'], [qubit_pair[1]])
        else:                             #same as above, but all operations are controlled
            if power%4 == 2:                #Z^2 
                qc.append(gates['cx'], [control_qubit, qubit_pair[0]])
            elif power%4 == 1:              #Z
                qc.append(gates['cx'], [control_qubit, qubit_pair[1]])
                qc.append(gates['cx'], [control_qubit, qubit_pair[0]])
                qc.append(gates['ccx'], [control_qubit, qubit_pair[1], qubit_pair[0]])
            elif power%4 == 3:              #Zdag
                qc.append(gates['ccx'], [control_qubit, qubit_pair[1], qubit_pair[0]])
                qc.append(gates['cx'], [control_qubit, qubit_pair[0]])
                qc.append(gates['cx'], [control_qubit, qubit_pair[1]])
        
    def X(self, qubit_pair, power = 1, control_qubit = -1, qc = None):
        if power%4 == 0:
            return
        if qc is None:
            qc = self.circuit
        if control_qubit == -1:        
            if power%4 == 2:                #X^2 
                qc.z(qubit_pair[1])
            elif power%4 == 1:              #X
                qc.z(qubit_pair[0])
                qc.s(qubit_pair[1])
            elif power%4 == 3:              #Xdag
                qc.z(qubit_pair[0])
                qc.sdg(qubit_pair[1])
        else:                             #same as above, but all operations are controlled
            if power%4 == 2:                #X^2 
                qc.cz(control_qubit, qubit_pair[1])
            elif power%4 == 1:              #X
                qc.cz(control_qubit, qubit_pair[0])
                qc.cp(pi/2, control_qubit, qubit_pair[1])
            elif power%4 == 3:              #Xdag
                qc.cz(control_qubit, qubit_pair[0])
                qc.cp(-pi/2, control_qubit, qubit_pair[1])
                
                
    def MagneticGS(self, qc = None, noisy = True):
        if qc is None:
            qc = self.circuit
        gates = self.gates(noisy)
        
        #all rows in parallel, three steps, sequential in the columns
        for x in range(self.Lx - 1):
            for i in range(3):
                #three iterations to perform in parallel these steps on each column but the last
                #i=0: hadamard the right link and act on right->bottom links
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
                        qc.append(gates['h'], [fourier_pair[0]])
                        qc.append(gates['h'], [fourier_pair[1]])
                        #qc.h(fourier_pair[0])
                        #qc.h(fourier_pair[1])
                        pair_index = 0
                    
                    pair = qubit_list[pair_index]
                    if pair_index > 1: qc.append(gates['x'], [pair[1]])      #X(4) if up or left link
                    qc.append(gates['ccx'], [pair[1], fourier_pair[1], pair[0]]) #CCX(4,2,3)
                    if pair_index > 1: qc.append(gates['x'], [pair[1]])      #X(4) if up or left link
                    qc.append(gates['cx'], [fourier_pair[0], pair[0]])      #CX(1,3)
                    qc.append(gates['cx'], [fourier_pair[1], pair[1]])      #CX(2,4)
                    """
                    if pair_index > 1: qc.x(pair[1])      #X(4) if up or left link
                    qc.ccx(pair[1], fourier_pair[1], pair[0]) #CCX(4,2,3)
                    if pair_index > 1: qc.x(pair[1])      #X(4) if up or left link
                    qc.cx(fourier_pair[0], pair[0])      #CX(1,3)
                    qc.cx(fourier_pair[1], pair[1])      #CX(2,4)
                    """
                                          
        #last column from top to bottom, last plaquette excluded       
        x = self.Lx - 1
        for y in range(self.Ly - 1, 0, -1):
            link_list = self.plaquette((x,y), from_zero = True)
            qubit_list = self.qubits_from_links(link_list)
            
            #fourier transform the lower link
            fourier_pair = qubit_list[0]
            qc.append(gates['h'], [fourier_pair[0]])
            qc.append(gates['h'], [fourier_pair[1]])
            
            for pair_index, pair in enumerate(qubit_list):
                if pair_index > 0:
                    if pair_index > 1: qc.append(gates['x'], [pair[1]])      #X(4) if up or left link
                    qc.append(gates['ccx'], [pair[1], fourier_pair[1], pair[0]]) #CCX(4,2,3)
                    if pair_index > 1: qc.append(gates['x'], [pair[1]])      #X(4) if up or left link
                    qc.append(gates['cx'], [fourier_pair[0], pair[0]])   #CX(1,3)
                    qc.append(gates['cx'], [fourier_pair[1], pair[1]])   #CX(2,4)
        return


class Z2(ToricCode):
    
    def __init__(self, size, pbc, N_ancillas):
        Lattice.__init__(self, size, pbc, spl = 2)
        self.N_ancillas = N_ancillas
        self.N_qubits = self.nlinks + self.N_ancillas
        self.circuit = QuantumCircuit(self.N_qubits, N_ancillas)
    
    def qubits_from_links(self, link_list):
        qubit_list = [link for link in link_list]
        return qubit_list
        
    def Z(self, qubit, power=1, control_qubit = -1, qc = None):
        if power%2==0:
            return
        if qc is None:
            qc = self.circuit
        if control_qubit == -1:        
            qc.x(qubit)
        else:                             
            qc.cx(control_qubit, qubit)
           
    def X(self, qubit, power=1, control_qubit = -1, qc = None):
        if power%2==0:
            return
        if qc is None:
            qc = self.circuit
        if control_qubit == -1:        
            qc.z(qubit)
        else:                             
            qc.cz(control_qubit, qubit)
                
                
    def MagneticGS(self, qc = None, noisy = True):
        if qc is None:
            qc = self.circuit
        h_label = cx_label = None
        if noisy is False:
            h_label = 'H_nl'
            cx_label = 'nl'
        h = ext.HGate(h_label)
        cx = ext.CXGate(cx_label)
            
        #all rows in parallel, three steps, sequential in the columns
        for x in range(self.Lx - 1):
            for i in range(3):
                #three iterations to perform in parallel these steps on each column but the last
                #i=0: hadamard the right link and act on right->bottom links
                #i=1: act on right->top links
                #i=2: act on right->left links (this is the only step that must
                #be done sequentially in the columns)
                for y in range(self.Ly):
                    link_list = self.plaquette((x,y), from_zero = True)
                    qubit_list = self.qubits_from_links(link_list)
                    fourier_qubit = qubit_list[1]
                    target_index = i+1
                    
                    #if first step:
                    if i==0:    
                        qc.append(h, [fourier_qubit])
                        #qc.h(fourier_qubit)                  
                        target_index = 0
                    
                    target_qubit = qubit_list[target_index]
                    qc.append(cx, [fourier_qubit, target_qubit])
                    #qc.cx(fourier_qubit, target_qubit)      #CX(1,2)
                                          
        #last column from top to bottom, last plaquette excluded       
        x = self.Lx - 1
        for y in range(self.Ly - 1, 0, -1):
            link_list = self.plaquette((x,y), from_zero = True)
            qubit_list = self.qubits_from_links(link_list)
            
            #fourier transform the lower link
            fourier_qubit = qubit_list[0]
            qc.append(h, [fourier_qubit])
            #qc.h(fourier_qubit)
            
            for target_index, target_qubit in enumerate(qubit_list):
                if target_index > 0:
                    qc.append(cx, [fourier_qubit, target_qubit])
                    #qc.cx(fourier_qubit, target_qubit)   #CX(1,2)