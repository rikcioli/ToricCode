# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 10:35:03 2022

@author: rikci
"""
import matplotlib.pyplot as plt
import json, codecs
import collections
import numpy as np
from IPython.display import display
import qiskit.providers.aer.noise as noise
from qiskit.visualization import plot_histogram

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
    prob_3 = 1 - ((1-prob_2)**13)*((1-prob_1)**12) # 3-qubit gate
    T1 = 180e3
    T2 = 180e3
    
    # Depolarizing quantum errors
    dep_gate = noise.depolarizing_error(prob_1, 1)
    dep_cx = noise.depolarizing_error(prob_2, 2)
    dep_ccx = noise.depolarizing_error(prob_3, 3)
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