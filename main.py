# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 18:31:08 2022

@author: rikci
"""

import ToricCode as tc
import numpy as np
from qiskit import transpile, IBMQ
from qiskit.providers.aer import QasmSimulator, AerSimulator
import qiskit.providers.aer.noise as noise
from qiskit.circuit.library import QFT
from qiskit.visualization import plot_histogram
from qiskit.tools import job_monitor
from IPython.display import display
import collections
import matplotlib.pyplot as plt
import json, codecs

def plot_counts(counts):
    #plot qiskit hist
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
    #prob_3 = 0.05
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
    return noise_model

def save_state(state, path = 'C:/Users/rikci/.spyder-py3/TESI/states/', name = 'Z4_3bit'):
    """
    state: np.array
    """
    lenght = len(state)
    dim = int(np.log2(lenght))
    re_sparse = [(index, re) for index, re in enumerate(state.real) if re!=0]
    im_sparse = [(index, im) for index, im in enumerate(state.imag) if im!=0]
    
    json.dump(re_sparse, codecs.open(path+name+'-real.json', 'w', encoding='utf-8'), 
          separators=(',', ':'), 
          sort_keys=True, 
          indent=4)
    json.dump(im_sparse, codecs.open(path+name+'-imag.json', 'w', encoding='utf-8'), 
          separators=(',', ':'), 
          sort_keys=True, 
          indent=4)
    json.dump(dim, codecs.open(path+name+'-dim.json', 'w', encoding='utf-8'),
          separators=(',', ':'), 
          sort_keys=True, 
          indent=4)
    return

def load_state(path = 'C:/Users/rikci/.spyder-py3/TESI/states/', name = 'Z4_3bit'):
    real_text = codecs.open(path+name+'-real.json', 'r', encoding='utf-8').read()
    imag_text = codecs.open(path+name+'-imag.json', 'r', encoding='utf-8').read()
    dim = int(codecs.open(path+name+'-dim.json', 'r', encoding='utf-8').read())
    re_sparse = json.loads(real_text)
    im_sparse = json.loads(imag_text)
    state = np.zeros(2**dim, complex)
    for index, re in re_sparse:
        state.real[index] = re
    for index, im in im_sparse:
        state.imag[index] = im
    return state



IBMQ.load_account()
provider = IBMQ.get_provider(hub='ibm-q')
backend = provider.get_backend('ibm_nairobi')
#noise_model = noise.NoiseModel.from_backend(backend)

# Initialize toric code and registers
N_reg = 3
test = tc.Z4((3,2), (True,True), N_reg)
pqc = test.circuit.copy()
#test.MagneticGS()
test.circuit.snapshot('final')

"""
ancilla = test.N_qubits-test.N_ancillas
reg_list = [ancilla+i for i in range(N_reg)]
test.circuit.h(reg_list)
"""

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
"""
test.init_dyons(low_charges = (2,1), up_charges = (3,1))
for i in range(N_reg):
    test.exchange_countclock(2**(i+1), control_qubit = reg_list[i])
"""

"""
test.circuit.append(QFT(N_reg, inverse=True), [test.circuit.qubits[reg] for reg in reg_list])
test.circuit.measure(reg_list, [i for i in range(N_reg)])
"""

# Draw the circuit and transpile it
qc_list = [test.circuit]
display(test.circuit.draw('mpl'))
simulator = AerSimulator(method = 'statevector')
tqc_list = transpile(qc_list, simulator)

"""
# Print circuit depth, size, width before and after transpile
print("\nTest circuit depth:", qc_list[0].depth())
print("Test circuit size:", qc_list[0].size())
print("Test circuit width:", qc_list[0].width())
print("\nTranspiled circuit depth:", tqc_list[0].depth())
print("Transpiled circuit size:", tqc_list[0].size())
print("Transpiled circuit width:", tqc_list[0].width())
"""

"""
# Run GS circuit on ideal simulator and save state on disk
N_shots = 1
job = simulator.run(tqc_list, job_name = "Ground State", shots = N_shots)
job_monitor(job)
print("\nRunning time {}s".format(job.result().time_taken))
snap1 = job.result().data()['snapshots']
state = snap1['statevector']['final'][0]
save_state(state, name = 'Z4_2bit')

#Run GS circuit on noisy simulator and measure average fidelity
N_noisy = 1
job_noisy = simulator.run(tqc_list, job_name = "Ground State noisy", shots = N_noisy)
job_monitor(job_noisy)
print("\nRunning time {}s".format(job.result().time_taken))
snap2 = job_noisy.result().data()['snapshots']
state_noisy_list = snap2['statevector']['final']
avg_fidelity = 0
for shot in state_noisy_list:
    fidelity = np.abs(np.dot(state.conjugate(), shot))**2
    avg_fidelity += fidelity/N_noisy
print("\nAverage fidelity with exact state is:", avg_fidelity)
"""

#Load ideal GS on a secondary circuit, and perform phase estimation
pqc.set_statevector(load_state(name = 'Z4_3bit'))
"""
ancilla = test.N_qubits-test.N_ancillas
reg_list = [ancilla+i for i in range(N_reg)]
pqc.h(reg_list)
#test.Z_string((1,1), (2,1), qc = pqc)
#test.X_string((1,0), (1,1), power = 1, qc = pqc)
#test.X_string((2,0), (2,1), power = 1, qc = pqc)
for i in range(N_reg):
    test.Bp((0,0), power=2*(2**i), qc = pqc, control_qubit=reg_list[i])
pqc.append(QFT(N_reg, inverse=True), [test.circuit.qubits[reg] for reg in reg_list])
pqc.measure(reg_list, [i for i in range(N_reg)])
"""
pqc.measure([18, 19, 20], [0, 1, 2])
tpqc = transpile(pqc, simulator)
display(tpqc.draw('mpl'))

N_shots_phase = 2000
job_phase = simulator.run(tpqc, job_name = "Phase Test", shots = N_shots_phase, noise_model = build_noise())
job_monitor(job_phase)
result_phase = job_phase.result()
print("\nRunning time {}s".format(result_phase.time_taken))
counts = result_phase.get_counts()
print("\nTotal counts are:", counts)
plot_counts(counts)
