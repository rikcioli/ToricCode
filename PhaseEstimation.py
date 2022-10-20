# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 18:31:08 2022

@author: rikci
"""

import ToricCodeNew as tc
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
import tools as ts

def plot_counts(counts):
    #plot qiskit hist
    counts = result.get_counts()
    display(plot_histogram(counts))
    #plot errorbar hist
    sorted_counts = dict(sorted(counts.items()))
    keylist = list(sorted_counts.keys())
    keycounts = list(sorted_counts.values())
    N_shots = sum(keycounts)
    yerr = np.sqrt(keycounts)/N_shots
    keycounts = np.array(keycounts)/N_shots
    plt.bar(keylist, keycounts, width = 0.5, yerr = yerr)
    plt.ylabel("Probabilities")
    plt.grid(axis = 'y', linestyle = '--')
    plt.xticks(rotation=70)
    plt.show()  

def show_results(result):
    #plot qiskit hist
    counts = result.get_counts()
    display(plot_histogram(counts))   
    # Extract memory for error estimates
    sorted_counts = dict(sorted(counts.items()))
    keylist = list(sorted_counts.keys())
    keycounts = list(sorted_counts.values())
    N_shots = sum(keycounts)
    memory = result.get_memory()
    N_bins = 10
    bin_width = N_shots//N_bins
    memory_binned = [memory[i*bin_width : (i+1)*bin_width] for i in range(N_bins)]
    hist_list = [dict(collections.Counter(memory_binned[i])) for i in range(N_bins)]
    
    dict_err = {}
    for key in keylist:
        for i in range(N_bins):
            if not key in hist_list[i]:
                hist_list[i][key] = 0
        key_results = [hist_list[i][key] for i in range(N_bins)]
        key_SDOM = np.std(key_results, ddof=1)/np.sqrt(N_bins)
        dict_err[key] = key_SDOM*N_bins
        
    yerr = list(dict_err.values())
    plt.bar(keylist, keycounts, width = 0.5, yerr = yerr)
    plt.show()

def build_noise():
    # Error probabilities
    prob_1 = 0.003  # 1-qubit gate
    prob_2 = 0.02   # 2-qubit gate
    prob_3 = 0.1
    #prob_3 = 1 - ((1-prob_2)**13)*((1-prob_1)**12) # 3-qubit gate
    
    # Depolarizing quantum errors
    error_1 = noise.depolarizing_error(prob_1, 1)
    error_2 = noise.depolarizing_error(prob_2, 2)
    error_3 = noise.depolarizing_error(prob_3, 3)
    
    # Add errors to noise model
    noise_model = noise.NoiseModel()
    noise_model.add_all_qubit_quantum_error(error_1, ['rz', 'sx', 'x', 'h', 'z', 's', 'sdg'])
    noise_model.add_all_qubit_quantum_error(error_2, ['cx', 'cp'])
    noise_model.add_all_qubit_quantum_error(error_3, ['ccx'])
    return noise_model

def build_noise_better():
    # Error probabilities
    prob_1 = 0.003  # 1-qubit gate
    prob_2 = 0.02   # 2-qubit gate
    prob_3 = 0.05
    #prob_3 = 1 - ((1-prob_2)**13)*((1-prob_1)**12) # 3-qubit gate
    T1 = 180e3
    T2 = 180e3
    
    # Depolarizing quantum errors
    dep_gate = noise.depolarizing_error(prob_1, 1)
    dep_cx = noise.depolarizing_error(prob_2, 2)
    dep_ccx = noise.depolarizing_error(prob_3, 3)
    # Thermal relaxation errors
    therm_gate = noise.thermal_relaxation_error(T1, T2, 50)   #estimated
    therm_cx = noise.thermal_relaxation_error(T1, T2, 536.9).expand(
             noise.thermal_relaxation_error(T1, T2, 536.9))     #exact
    therm_reset = noise.thermal_relaxation_error(T1, T2, 4000)    #estimated
    therm_readout = noise.thermal_relaxation_error(T1, T2, 4924.4)    #exact
    
    # Add errors to noise model
    noise_model = noise.NoiseModel()
    noise_model.add_all_qubit_quantum_error(dep_gate.compose(therm_gate), ['id', 'rz', 'sx', 'x', 'h', 'z', 's', 'sdg'])
    noise_model.add_all_qubit_quantum_error(dep_cx.compose(therm_cx), ['cx', 'cp'])
    noise_model.add_all_qubit_quantum_error(dep_ccx, ['ccx'])
    noise_model.add_all_qubit_quantum_error(therm_reset, ['reset'])
    noise_model.add_all_qubit_quantum_error(therm_readout, ['measure'])
    return noise_model

#IBMQ.enable_account('267761afd846893dec77bf06dc487d6c7b9569ed20605f88fbddf79e6fb4d149dcfb04da0e4c78994302537adbb9da62c7cfc1d5341f9d4a4c8b416a953bbc9b')
IBMQ.load_account()
provider = IBMQ.get_provider(hub='ibm-q')
backend = provider.get_backend('ibm_nairobi')

#provider = IBMQ.get_provider(hub='ibm-q-cern')
#backend = provider.get_backend('ibm_cairo')
noise_model = noise.NoiseModel.from_backend(backend)
simulator = AerSimulator(method = "statevector")

# Initialize toric code and registers
N_reg = 3
test = tc.Z3((3,2), (True, True), N_reg)
test.MagneticGS()
ancilla = test.N_qubits-test.N_ancillas
reg_list = [ancilla+i for i in range(N_reg)]
test.circuit.h(reg_list)

#Initialize e and m particles and rotate them

test.Z_string((1,1), (2,1))
test.X_string((1,0), (1,1), power = 1)
#test.X_string((2,0), (2,1), power = 1)
for i in range(N_reg):
    test.Bp((0,0), power=2*(2**i), control_qubit=reg_list[i])

#Access non trivial sector and measure it with 't Hooft loop
"""
test.Z_string((0,1), (3,1), power = 2)
test.X_string((1,0), (2,0), control_qubit=ancilla)
test.X_string((1,1), (2,1), control_qubit=ancilla)
"""
#Initialize dyons and exchange them an arbitrary number of times
"""
test.init_dyons(low_charges = (2,1), up_charges = (3,1))
for i in range(N_reg):
    test.exchange_countclock(2**(i+1), control_qubit = reg_list[i])
"""

test.circuit.append(QFT(N_reg, inverse=True), [test.circuit.qubits[reg] for reg in reg_list])
test.circuit.measure(reg_list, [i for i in range(N_reg)])

display(test.circuit.draw('mpl'))

# Transpile the circuit multiple times for stochastic swap
qc_list = [test.circuit.copy() for i in range(10)]
tqc_list = transpile(qc_list, 
                     simulator, 
                     optimization_level=1, 
                     )

# Find best circuit
depths = [tqc.depth() for tqc in tqc_list]
min_depth_index = min(range(len(depths)), key=depths.__getitem__)
tqc = tqc_list[min_depth_index]
display(tqc.draw('mpl'))


# Print circuit depth, size, width before and after transpile
print("\nTest circuit depth:", test.circuit.depth())
print("Test circuit size:", test.circuit.size())
print("Test circuit width:", test.circuit.width())
print("\nTranspiled circuit depth:", tqc.depth())
print("Transpiled circuit size:", tqc.size())
print("Transpiled circuit width:", tqc.width())

# Run circuit on simulator/backend and get results
N_shots = 2000
job = simulator.run(tqc, job_name = "Toric Test", shots = N_shots)
job_monitor(job)
result = job.result()

# Extract counts
counts = result.get_counts()
print("\nTotal counts are:", counts)
ts.plot_counts(counts)

# Print time taken
print("\nRunning time {}s".format(result.time_taken))


