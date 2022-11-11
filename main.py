# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 18:31:08 2022

@author: rikci
"""

import math
import ToricCodeNewFinal as tc
import tools as ts
import scipy.sparse as sp
from qiskit import transpile, IBMQ, QuantumCircuit
from qiskit.providers.aer import AerSimulator
import qiskit.providers.aer.noise as noise
from qiskit.circuit.library import QFT
from qiskit.visualization import plot_histogram
from qiskit.tools import job_monitor
from IPython.display import display


def run_heavy(tqc, shots):
    job = simulator.run(tqc, job_name = "GS Check", shots = shots)
    job_monitor(job)
    result = job.result()
    print("Running time {}s\n".format(result.time_taken))
    display(plot_histogram(result.get_counts()))
        

def check_GS(coo_state, simulator, N_reg = 3, plaq = False):
    """Check GS copy by evaluating all plaquette and star eigenvalues"""
    
    state = coo_state.toarray()[0]
    dim = int(math.log2(len(state)))
    qc = QuantumCircuit(dim, dim)
    qc.set_statevector(state)
    
    ancilla = dim - N_reg
    reg_list = [ancilla+i for i in range(N_reg)]
    qc.reset(reg_list)
    qc.h(reg_list)
    
    if plaq is True:
        for i in range(N_reg):
            test.Bp((0,0), power=3*(2**i), qc = qc, control_qubit=reg_list[i])
        qc.append(QFT(N_reg, inverse=True), [qc.qubits[reg] for reg in reg_list])
        qc.measure(reg_list, [j for j in range(N_reg)])
        qc_list = [qc]
    
    else:     
        qc_list = [qc.copy() for i in range(6)]
        
        for qc_index, GS in enumerate(qc_list):
            x = 0
            y = 1
            if qc_index%3 == 2: x = 1
            if qc_index%3 == 0: y = 0
            for i in range(N_reg):
                if qc_index < 3:
                    test.Bp((x,y), power=3*(2**i), qc = GS, control_qubit=reg_list[i])
                else:
                    test.As((x,y), power=3*(2**i), qc = GS, control_qubit=reg_list[i])
            GS.append(QFT(N_reg, inverse=True), [qc.qubits[reg] for reg in reg_list])
            GS.measure(reg_list, [j for j in range(N_reg)])
    
    N_shots_check = 10
    tqc_list = transpile(qc_list, simulator)
    for tGS in tqc_list:
        run_heavy(tqc = tGS, shots = N_shots_check)
        
        
def runSave(tqc, final_perm = None, filename = None, **kwargs):
    job = simulator.run(tqc, **kwargs)
    job_monitor(job)
    print("Running time {}s\n".format(job.result().time_taken))
    state = job.result().data()['snapshots']['statevector']['final'][0]
    state = sp.coo_matrix(ts.cleanState(state))
    if final_perm is not None:
        state = ts.swapSparseState(state, final_perm)
    if filename is not None:
        sp.save_npz(filename, state)
    return state


#IBMQ.enable_account('267761afd846893dec77bf06dc487d6c7b9569ed20605f88fbddf79e6fb4d149dcfb04da0e4c78994302537adbb9da62c7cfc1d5341f9d4a4c8b416a953bbc9b')
#IBMQ.load_account()
provider = IBMQ.get_provider(hub='ibm-q-cern')
backend = provider.get_backend('ibm_cairo')
noise_model = noise.NoiseModel.from_backend(backend)
simulator = AerSimulator(method = 'statevector')

# Initialize toric code and registers s.t. number of qubits matches with backend
test = tc.Z3((3,2), (False, False), backend = backend)
test.MagneticGS_2()
display(test.circuit.draw('mpl'))

"""
layout2x2 = [4, 5, 6, 7, 14, 15, 12, 13, 0, 1, 2, 3, 10, 11, 8, 9] # use with MagneticGS
[layout2x2.append(i) for i in range(16, 27)]
layout2x2bis = [8, 9, 10, 11, 2, 3, 0, 1, 4, 5, 6, 7, 14, 15, 12, 13] # use with MagneticGS_2
[layout2x2bis.append(i) for i in range(16, 27)]
layout2plaq = [4, 5, 10, 11, 2, 3, 8, 9, 0, 1, 6, 7, 12, 13] # use with 2 plaq, MagneticGS_2
[layout2plaq.append(i) for i in range(14, 27)]

qc = transpile(test.circuit,
               initial_layout = layout2x2bis)
display(qc.draw('mpl'))
"""

# Measure all qubits to identify final permutations. The measurements will be removed later
test.circuit.measure([i for i in range(test.N_qubits)], [i for i in range(test.N_qubits)])

# Transpile the circuit multiple times for stochastic swap
qc_list = [test.circuit.copy() for i in range(50)]
tqc_list = transpile(qc_list, 
                     backend, 
                     optimization_level=3,
                     )

# Find best circuit
depths = [tqc.depth() for tqc in tqc_list]
min_depth_index = min(range(len(depths)), key=depths.__getitem__)
tqc = tqc_list[min_depth_index]
display(tqc.draw('mpl'))

# Extract virtual to real global permutation
final_perm = [-1]*27
for op, qubits, clbits in tqc.data:
       if op.name == 'measure':
           final_perm[clbits[0].index] = qubits[0].index

# Remove final measurements and barrier
for i in range(test.N_qubits+1):
    tqc.data.pop(len(tqc.data)-1)

# Save snapshot and print final circuit
tqc.snapshot('final')
display(tqc.draw('mpl'))
print("\nTest circuit depth:", test.circuit.depth())
print("Test circuit size:", test.circuit.size())
print("Test circuit width:", test.circuit.width())
print("\nTranspiled circuit depth:", tqc.depth())
print("Transpiled circuit size:", tqc.size())
print("Transpiled circuit width:", tqc.width())

# Run GS circuit on ideal simulator, extract statevector as sparse array and save it
state_ideal = runSave(tqc,
                      final_perm,
                      filename = 'C:/Users/rikci/.spyder-py3/TESI/SciPy States/test/ideal',
                      shots = 1)

#Run GS circuit copy on noisy simulator and measure average fidelity
N_states_noisy = 1000
states = [runSave(tqc, 
                  final_perm, 
                  filename = 'C:/Users/rikci/.spyder-py3/TESI/SciPy States/test/'+str(i+1), 
                  shots = 1,
                  #noise_model = ts.build_noise()
                  ) for i in range(N_states_noisy)]

fidelity = ts.evaluateFidelity(states, state_ideal)

print("Fidelity with exact state is:", fidelity)
