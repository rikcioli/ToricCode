# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 12:11:05 2022

@author: rikci
"""

import numpy as np
from qiskit import QuantumCircuit, transpile, IBMQ
from qiskit.providers.aer import QasmSimulator
from qiskit.visualization import plot_histogram
import qiskit.providers.aer.noise as noise
from qiskit.tools import job_monitor
from qiskit.circuit.library import QFT
from IPython.display import display
import phases

#IBMQ.load_account() # Load account from disk
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


simulator = QasmSimulator()

N_shots = 20000
N_registers = 2
phase = 1/4

qc = QuantumCircuit(N_registers+1, N_registers+1)
reg_list = [i for i in range(1, N_registers+1)]
qc.h(reg_list)
qc.x(0)
for qubit_index in reg_list:
    for j in range(2**(qubit_index-1)):
        qc.cp(2*np.pi*phase, qubit_index, 0)
        
qc.append(QFT(N_registers, inverse = True), [qc.qubits[i] for i in reg_list])
qc.measure_all()

"""
qc = QuantumCircuit(2, 2)
theta = -2*np.arccos(1/np.sqrt(3))
qc.ry(theta, 0)
qc.ch(0, 1)
qc.x(0)
qc.measure([0,1], [0,1])
"""

tqc = transpile(qc, backend)
display(qc.draw())
display(tqc.draw())

job = simulator.run(tqc, job_name = "simple_test", shots=N_shots, memory=True, noise_model = noise_model)
job_monitor(job)
results = job.result()
counts = results.get_counts()

print("\nTotal counts are:", counts)
print("\nTest circuit depth:", qc.depth())
print("Test circuit size:", qc.size())
print("Test circuit width:", qc.width())
print("\nTranspiled circuit depth:", tqc.depth())
print("Transpiled circuit size:", tqc.size())
print("Transpiled circuit width:", tqc.width())

display(plot_histogram(counts))