# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 17:29:04 2022

@author: rikci
"""

import numpy as np

def evaluate_phase(counts):   
    #Given the counts of a measurement of the ancilla, returns the relative phase
    #between the 2 Z eigenstates   
    phase = 0
    keys = list(counts.keys())   
    if len(keys) == 1:
        if int(keys[0]) > 0:
            phase = np.pi
    else:
        extracted_values = list(counts.values())
        if int(keys[0]) == 0:
            values = extracted_values
        else:
            values = [extracted_values[1], extracted_values[0]]
        cos_phi = (values[0] - values[1])/(values[0] + values[1])
        phase = np.arccos(cos_phi)
    return phase

def evaluate_sin(counts):
    sin_phi = 1
    keys = list(counts.keys())
    if len(keys) == 1:
        if int(keys[0]) > 0:
            sin_phi = -1
    else:
        extracted_values = list(counts.values())
        if int(keys[0]) == 0:
            values = extracted_values
        else:
            values = [extracted_values[1], extracted_values[0]]
        sin_phi = (values[0] - values[1])/(values[0] + values[1])
    return sin_phi

def binning(memory, N_shots, N_bins):
    """
    Organize the full data into bins, each bin simulating a single experiment with
    number of shots: N_shots/N_bins

    Parameters
    ----------
    memory : TYPE
        DESCRIPTION.
    N_shots : TYPE
        DESCRIPTION.
    N_bins : TYPE
        DESCRIPTION.

    Returns
    -------
    binned_data : Dictionary
        Key: bin number. Value: Dictionary containing the outcomes of the 
        experiment associated to the given bin number.

    """  

    binned_data = {}
    bin_width = N_shots//N_bins
    for i in range(N_bins):
        N_0 = 0
        N_1 = 0
        for j in range(bin_width):
            index = i*bin_width + j
            if int(memory[index]) == 0:
                N_0 += 1
            else:
                N_1 += 1
        binned_data[i] = {}
        binned_data[i]["0"] = N_0
        binned_data[i]["1"] = N_1
    return binned_data

def avg_phase(binned_data):
    """
    

    Parameters
    ----------
    binned_data : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    #Given a dictionary of experiments, extract a value of theta from each experiment
    thetas = []
    results = list(binned_data.values())  #extract list of experiments from dict of experiments
    for single_exp in results:
        thetas.append(evaluate_phase(single_exp))
    thetas = np.array(thetas)
    mean = thetas.mean()
    SDOM = np.std(thetas, ddof=1)/np.sqrt(len(thetas))
    return (mean, SDOM)

def avg_sin(binned_data):
    #Given a dictionary of experiments, extract a value of sin(theta) from each experiment
    sin_thetas = []
    results = list(binned_data.values())  #convert dict of experiments to a list of experiments
    for single_exp in results:
        sin_thetas.append(evaluate_sin(single_exp))
    sin_thetas = np.array(sin_thetas)
    mean = sin_thetas.mean()
    SDOM = np.std(sin_thetas, ddof=1)/np.sqrt(len(sin_thetas))
    return (mean, SDOM)

def results(job_result, N_shots, N_bins = 50, N_cos_circuits = 1, N_sin_circuits = 1):
    N_circuits = N_cos_circuits + N_sin_circuits
    counts_list = [job_result.get_counts(i) for i in range(N_circuits)]
    memory_list = [job_result.get_memory(i) for i in range(N_circuits)]
    data_list = [binning(memory, N_shots = N_shots, N_bins = N_bins) for memory in memory_list]
    phase_list = []
    for i, binned_data in enumerate(data_list):
        if i < N_cos_circuits:
            print("\n\nCIRCUIT N:", i+1, "- COS")
            phase_list.append(avg_phase(binned_data))
            phase_result = phase_list[i]
            print("\nTotal counts are:", counts_list[i])
            print("Total phase is:", evaluate_phase(counts_list[i]))
            print("\nAverage phase is:", phase_result[0])
            print("Error on phase is:", phase_result[1])
            if phase_result[0]!=0: 
                print("Relative uncertainty:", 100*phase_result[1]/phase_result[0], "%")
        else:
            print("\n\nCIRCUIT N:", i - N_cos_circuits + 1, "- SIN")
            phase_list.append(avg_sin(binned_data))
            sin_result = phase_list[i]
            print("\nTotal counts are:", counts_list[i])
            print("Total sin is:", evaluate_sin(counts_list[i]))
            print("\nAverage sin is:", sin_result[0])
            print("Error on sin is:", sin_result[1])
            if sin_result[0]!=0: 
                print("Relative uncertainty:", 100*sin_result[1]/np.abs(sin_result[0]), "%")
    return phase_list

