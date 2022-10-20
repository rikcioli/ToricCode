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
    QUBIT 0 = 2*LINK
    QUBIT 1 = 2*LINK + 1
Except for the registers, which are the qubits ranging from N_qubits-N_ancillas to N_qubits
"""


class ToricCode(Lattice):
    
    def __init__(self, size, pbc, spl, N_ancillas = 0, backend = None):
        super().__init__(size, pbc, spl)
        self.N_ancillas = N_ancillas
        """
        self.N_qubits = 2*self.nlinks + self.N_ancillas
        self.circuit = QuantumCircuit(self.N_qubits, self.N_qubits)
        self.reg_list = [(2*self.nlinks + i) for i in range(N_ancillas)]
        """
    def _qubits_from_link(self, link):
        """Returns the qubits associated to a given link"""
        qubits = (2*link, 2*link+1)
        return qubits
    
    def _gates(self, noisy = True):
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
        

    def _Z(self, link, power = 1, control_qubit = -1, qc = None, noisy = True):
        return
        
    def _X(self, link, power = 1, control_qubit = -1, qc = None, noisy = True):
        return         
    
    def _elecToMagnGS(self, link, qc = None, noisy = True):
        return
    
    def _CU(self, control_link, target_link, inverse = False, qc = None, noisy = True):
        return
            
    def MagneticGS(self, qc = None, noisy = True):
        if self.nlinks_y > self.nlinks_x:
            raise RuntimeError("nlinks_y is greater than nlinks_x. To correctly initialize the Magnetic GS, please use a lattice with nlinks_x >= nlinks_y")
            
        #all rows in parallel, three steps, sequential in the columns
        for x in range(self.Lx - 1):
            for i in range(3):
                #three iterations to perform in parallel these steps on each column but the last
                #i=0: fourier transform the right link and act on right->bottom links
                #i=1: act on right->top links
                #i=2: act on right->left links (this is the only step that must
                #be done sequentially in the columns)
                inverse = False
                for y in range(self.nlinks_y):                    
                    link_list = self.plaquette((x,y), from_zero = True)
                    control_link = link_list[1] #right link
                    target_index = i+1
                    
                    #if first step:
                    if i==0:    
                        self._elecToMagnGS(control_link, qc, noisy)
                        target_index = 0
                    
                    target_link = link_list[target_index]
                    if target_index > 1: inverse = True
                    self._CU(control_link, target_link, inverse, qc, noisy)
                                          
        if self.pbc_x == True:
            #last column from top to bottom, last plaquette excluded       
            x = self.Lx - 1
            for y in range(self.Ly - 1, 0, -1):
                link_list = self.plaquette((x,y), from_zero = True)
                control_link = link_list[0]
                
                #fourier transform the lower link
                self._elecToMagnGS(control_link, qc, noisy)
    
                for target_index, target_link in enumerate(link_list):
                    inverse = False
                    if target_index > 0:
                        if target_index > 1:
                            inverse = True
                        self._CU(control_link, target_link, inverse, qc, noisy)
        return
    
                    
    def Z_string(self, site0, site1, power = 1, **kwargs):        
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
        for link in path:
            self._Z(link, power, **kwargs) 
                    
    def X_string(self, site0, site1, power = 1, **kwargs):
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
        for link in path:
            self._X(link, power, **kwargs)
     
    def Bp(self, site, power = 1, **kwargs):
        if power%self.spl == 0:
            return
        x, y = site
        sites = [(x, y), (x+1, y), (x+1, y+1), (x, y+1)]
        for i in range(4):
            self.Z_string(sites[i], sites[(i+1)%4], power, **kwargs)
    
    def As(self, site, power = 1, **kwargs):
        if power%self.spl == 0:
            return
        x, y = site
        endsites = [(x+1, y), (x, y+1), (x-1, y), (x, y-1)]
        for endsite in endsites:
            self.X_string(site, endsite, power, **kwargs)
        
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
    
    def __init__(self, size, pbc, N_ancillas = 0, backend = None):
        super().__init__(size, pbc, 3, N_ancillas, backend)
        N_qubits = 2*self.nlinks + self.N_ancillas
        if backend is None:
            self.N_qubits = N_qubits
        else:
            N_qubits_backend = backend.configuration().n_qubits
            if N_qubits_backend < N_qubits:
                raise RuntimeError('Number of qubits required exceeds available on given backend')
            else:
                self.N_qubits = N_qubits_backend
        self.circuit = QuantumCircuit(self.N_qubits, self.N_qubits)
        self.reg_list = [(2*self.nlinks + i) for i in range(N_ancillas)]
        
    
    def _Z(self, link, power = 1, control_qubit = -1, qc = None, noisy = True):
        if power%3 == 0:
            return
        if qc is None:
            qc = self.circuit
        gates = self._gates(noisy)
        qubit_pair = self._qubits_from_link(link)
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
    
    def _X(self, link, power = 1, control_qubit = -1, qc = None, noisy = True):
        #NOISELESS OPTION HAS YET TO BE IMPLEMENTED
        if power%3 == 0:
            return
        if qc is None:
            qc = self.circuit
        qubit_pair = self._qubits_from_link(link)
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
                
    def _elecToMagnGS(self, link, qc = None, noisy = True):
        """Convert link from electric GS to magnetic GS"""
        if qc is None:
            qc = self.circuit
        gates = self._gates(noisy)
        qubit_pair = self._qubits_from_link(link)
        qc.append(gates['ry'], [qubit_pair[0]])
        qc.append(gates['ch'], [qubit_pair[0], qubit_pair[1]])
        qc.append(gates['x'], [qubit_pair[0]])
        
    def _CU(self, control_link, target_link, inverse, qc = None, noisy = True):
        """Z3 entanglement operation between two links"""
        if qc is None:
            qc = self.circuit
        gates = self._gates(noisy)
        control_qubits = self._qubits_from_link(control_link)
        
        qc.append(gates['cx'], [control_qubits[0], control_qubits[1]])
        if inverse == False:
            self._Z(target_link, power = 1, control_qubit = control_qubits[0], noisy=noisy)
            self._Z(target_link, power = 1, control_qubit = control_qubits[1], noisy=noisy)
        else:
            self._Z(target_link, power = 2, control_qubit = control_qubits[1], noisy=noisy)
            self._Z(target_link, power = 2, control_qubit = control_qubits[0], noisy=noisy)
        qc.append(gates['cx'], [control_qubits[0], control_qubits[1]])
       
                
                
class Z4(ToricCode):
    
    def __init__(self, size, pbc, N_ancillas = 0, backend = None):
        super().__init__(size, pbc, 4, N_ancillas, backend)
        N_qubits = 2*self.nlinks + self.N_ancillas
        if backend is None:
            self.N_qubits = N_qubits
        else:
            N_qubits_backend = backend.configuration().n_qubits
            if N_qubits_backend < N_qubits:
                raise RuntimeError('Number of qubits required exceeds available on given backend')
            else:
                self.N_qubits = N_qubits_backend
        self.circuit = QuantumCircuit(self.N_qubits, self.N_qubits)
        self.reg_list = [(2*self.nlinks + i) for i in range(N_ancillas)]
        
        
    def _Z(self, link, power = 1, control_qubit = -1, qc = None, noisy = True):
        if power%4 == 0:
            return
        if qc is None:
            qc = self.circuit
        gates = self._gates(noisy)
        qubit_pair = self._qubits_from_link(link)
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
        
    def _X(self, link, power = 1, control_qubit = -1, qc = None, noisy = True):
        #NOISELESS OPTION HAS YET TO BE IMPLEMENTED
        if power%4 == 0:
            return
        if qc is None:
            qc = self.circuit
        qubit_pair = self._qubits_from_link(link)
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
                
    def _elecToMagnGS(self, link, qc = None, noisy = True):
        """Convert link from electric GS to magnetic GS"""
        if qc is None:
            qc = self.circuit
        gates = self._gates(noisy)
        qubit_pair = self._qubits_from_link(link)
        qc.append(gates['h'], [qubit_pair[0]])
        qc.append(gates['h'], [qubit_pair[1]])
        
        
    def _CU(self, control_link, target_link, inverse, qc = None, noisy = True):
        """Z4 entanglement operation between two links"""
        if qc is None:
            qc = self.circuit
        gates = self._gates(noisy)
        control_qubits = self._qubits_from_link(control_link)
        target_qubits = self._qubits_from_link(target_link)
        
        # X(4) if up or left link
        if inverse == True: qc.append(gates['x'], [target_qubits[1]]) 
        
        # CCX(4,2,3)
        qc.append(gates['ccx'], [target_qubits[1], control_qubits[1], target_qubits[0]]) 
        
        # X(4) if up or left link
        if inverse == True: qc.append(gates['x'], [target_qubits[1]]) 
        
        # CX(1,3) and #CX(2,4)       
        qc.append(gates['cx'], [control_qubits[0], target_qubits[0]])      
        qc.append(gates['cx'], [control_qubits[1], target_qubits[1]])      


class Z2(ToricCode):
    
    def __init__(self, size, pbc, N_ancillas = 0, backend = None):
        super().__init__(size, pbc, 2, N_ancillas, backend)
        N_qubits = self.nlinks + self.N_ancillas
        if backend is None:
            self.N_qubits = N_qubits
        else:
            N_qubits_backend = backend.configuration().n_qubits
            if N_qubits_backend < N_qubits:
                raise RuntimeError('Number of qubits required exceeds available on given backend')
            else:
                self.N_qubits = N_qubits_backend
        self.circuit = QuantumCircuit(self.N_qubits, self.N_qubits)
        self.reg_list = [(self.nlinks + i) for i in range(N_ancillas)]
        
        
    def _Z(self, link, power=1, control_qubit = -1, qc = None, noisy = True):
        if power%2==0:
            return
        if qc is None:
            qc = self.circuit
        if control_qubit == -1:        
            qc.x(link)
        else:                             
            qc.cx(control_qubit, link)
           
    def _X(self, link, power=1, control_qubit = -1, qc = None, noisy = True):
        if power%2==0:
            return
        if qc is None:
            qc = self.circuit
        if control_qubit == -1:        
            qc.z(link)
        else:                             
            qc.cz(control_qubit, link)
            
    def _elecToMagnGS(self, link, qc = None, noisy = True):
        """Convert link from electric GS to magnetic GS"""
        if qc is None:
            qc = self.circuit
        gates = self._gates(noisy)
        qc.append(gates['h'], [link])
    
    def _CU(self, control_link, target_link, inverse, qc = None, noisy = True):
        """Z4 entanglement operation between two links"""
        if qc is None:
            qc = self.circuit
        gates = self._gates(noisy)
        # CX(1,2)
        qc.append(gates['cx'], [control_link, target_link])
                