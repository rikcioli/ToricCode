# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 18:31:08 2022

@author: rikci
"""

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.providers.aer import QasmSimulator
from IPython.display import display

N_links = 10
N_ancillas = 1
N_qubits = N_links+N_ancillas
N_shots = 10000
circuit = QuantumCircuit(N_qubits,N_qubits)

def ElectricGS():
    for i in range(N_links):
        circuit.h(i)

def MagneticGS():
    """
    Initializes magnetic ground state, by applying (1+Bp) on each plaquette.
    The application proceeds from left to right, using the right link as target of
    the other 3. All the non-target qubits need to be initialized to eigenstates
    of X with eigenvalue +1.    
    """
    target = 3
    for i in range(N_links-1):       
        if i%3 == 0:
            target = (i+3)
        else:
            circuit.h(i)
        circuit.cx(i, target)
    
def Z_string(links):
    for i in links:
        circuit.z(i)
    
def X_string(links):
    for i in links:
        circuit.x(i)
        
def evaluate_phase(counts):
    """
    Given the counts of a measurement of the ancilla, returns the relative phase
    between the 2 Z eigenstates
    """
    values = list(counts.values())
    keys = list(counts.keys())
    phase = 0
    if len(keys) == 1:
        if int(keys[0]) > 0:
            phase = np.pi
    else:
        cos_phi = (values[0] - values[1])/N_shots
        phase = np.arccos(cos_phi)
    return phase


#Prepare magnetic GS
MagneticGS()

#Initialize ancilla
ancilla = N_qubits-1
circuit.h(ancilla)

#Initialize Z string and X string
Z_string([5])
X_string([3,6])

#Rotation of e around m controlled by the ancilla
circuit.cz(ancilla, 2)
circuit.cz(ancilla, 0)
circuit.cz(ancilla, 1)
circuit.cz(ancilla, 3)

#Hadamard back the ancilla and measure it
circuit.h(ancilla)
circuit.measure(N_qubits-1, N_qubits-1)

#Draw the circuit and run it through the simulator with a given number of shots
display(circuit.draw())
simulator = QasmSimulator()
compiled_circuit = transpile(circuit, simulator)
job = simulator.run(compiled_circuit, shots=N_shots)
result = job.result()
counts = result.get_counts(compiled_circuit)
print("\nTotal counts are:",counts)

#Extract relative phase from the counts
phase = evaluate_phase(counts)
print("\nTotal phase is:", phase)
