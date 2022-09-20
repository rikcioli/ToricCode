# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 18:31:08 2022

@author: rikci
"""

import ToricCode as tc
import numpy as np
from qiskit import transpile, IBMQ
from qiskit.providers.aer import AerSimulator
import qiskit.providers.aer.noise as noise
from qiskit.circuit.library import QFT
from qiskit.visualization import plot_histogram
from qiskit.tools import job_monitor
from IPython.display import display
import collections
import matplotlib.pyplot as plt

IBMQ.load_account()
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


# Initialize Z3 toric code and registers
N_reg = 2
test = tc.Z4((3,2), (True,True), N_reg)
test.MagneticGS()
ancilla = test.N_qubits-test.N_ancillas
reg_list = [ancilla+i for i in range(N_reg)]
test.circuit.h(reg_list)

#Initialize e and m particles and rotate them
"""
test.Z_string((1,1), (2,1))
test.X_string((1,0), (1,1), power = 1)
test.X_string((2,0), (2,1), power = 1)
for i in range(N_reg):
    test.Bp((0,0), power=3*(2**i), control_qubit=reg_list[i])
"""
#Access non trivial sector and measure it with 't Hooft loop
"""
test.Z_string((0,1), (3,1), power = 2)
test.X_string((1,0), (2,0), control_qubit=ancilla)
test.X_string((1,1), (2,1), control_qubit=ancilla)
"""
#Initialize dyons and exchange them an arbitrary number of times

test.init_dyons(low_charges = (2,1), up_charges = (3,1))
for i in range(N_reg):
    test.exchange_countclock(2**(i+1), control_qubit = reg_list[i])



#test.circuit.h(reg_list)
test.circuit.append(QFT(N_reg, inverse=True), [test.circuit.qubits[reg] for reg in reg_list])
test.circuit.measure(reg_list, [i for i in range(N_reg)])
qc_list = [test.circuit]


# Draw the circuit and transpile it
display(test.circuit.draw('mpl'))
simulator = AerSimulator(method = "statevector")
tqc_list = transpile(qc_list, simulator)

# Print circuit depth, size, width before and after transpile
print("\nTest circuit depth:", qc_list[0].depth())
print("Test circuit size:", qc_list[0].size())
print("Test circuit width:", qc_list[0].width())
print("\nTranspiled circuit depth:", tqc_list[0].depth())
print("Transpiled circuit size:", tqc_list[0].size())
print("Transpiled circuit width:", tqc_list[0].width())

# Run circuit on simulator/backend and get results
N_shots = 2000
job = simulator.run(tqc_list, job_name = "Toric Test", shots = N_shots, memory = True)
job_monitor(job)
result = job.result()

# Extract counts
counts = result.get_counts()
print("\nTotal counts are:", counts)
display(plot_histogram(counts))

# Extract memory for error estimates
memory = result.get_memory()
N_bins = 10
bin_width = N_shots//N_bins
memory_binned = [memory[i*bin_width : (i+1)*bin_width] for i in range(N_bins)]
hist_list = [dict(collections.Counter(memory_binned[i])) for i in range(N_bins)]

sorted_counts = dict(sorted(counts.items()))
keylist = list(sorted_counts.keys())

dict_err = {}
for key in keylist:
    for i in range(N_bins):
        if not key in hist_list[i]:
            hist_list[i][key] = 0
    key_results = [hist_list[i][key] for i in range(N_bins)]
    key_SDOM = np.std(key_results, ddof=1)/np.sqrt(N_bins)
    dict_err[key] = key_SDOM*N_bins
    

plt.bar(keylist, list(sorted_counts.values()), width = 0.5, yerr = list(dict_err.values()))
plt.show()

# Print time taken
print("\nRunning time {}s".format(result.time_taken))


    

