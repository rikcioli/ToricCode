# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 12:11:05 2022

@author: rikci
"""

import numpy as np
import scipy.sparse as sp
from qiskit import QuantumCircuit, transpile, IBMQ
from qiskit.providers.aer import QasmSimulator, AerSimulator
from qiskit.visualization import plot_histogram, plot_state_city, plot_state_hinton
from qiskit.providers.aer.noise import NoiseModel
from qiskit.tools import job_monitor
from IPython.display import display
import phases
import tools as ts

def swapBits(x, p1, p2):
    """Swap bits of positions p1 and p2 in number x"""
    xor = (((x >> p1) ^ (x >> p2)) & 1)
    return x ^ ( (xor << p1) | (xor << p2))

def findPerm(target_list):
    """target_list[i] is a value associated to index i. The index represents the
    starting configuration. The value represents the final configuration after
    the permutation. So for example, target_list[0] = 3 means the object of 
    position 0 will become the object of position 3. 
    This is opposite to qiskit initial_layout assignments, where the starting qubits
    take the place of the end qubits."""
    permutations = []
    graph = {
        start: end
        for start, end in enumerate(target_list)
        }
    for start, end in graph.items():
        check_presence = [start in perm for perm in permutations]
        if True not in check_presence:
            if end != start:
                perm = [start, end]      
                newstep = graph[end]
                while newstep != start:
                    perm.append(newstep)
                    newstep = graph[newstep]
                permutations.append(perm.copy())
    return permutations

def permToSwap(permutations_list):
    swaps = []
    for perm in permutations_list:
        for i in range(len(perm)-1):
            swaps.append(perm[i:i+2])
    return swaps
            

def swapState(state, assign_list):
    """Swap qubits in a given state (state is an array of complex). The assign
    list works this way: the index runs over the starting positions, while the values
    represent the end configuration, after the permutation. assign_list[i] = j
    means that i COPIES j. It DOES NOT take its place, which is how Qiskit instead
    operates. If you want to recover the original Qiskit statevector as it was
    before the transpilation, simply put the initial_layout of the transpiler into
    assign_list."""
    permutations = findPerm(assign_list)
    swaps = permToSwap(permutations)
    state_sparse = sy.coo_matrix(state)
    for index in range(state_sparse.nnz): 
        for swap in swaps:
            state_sparse.col[index] = swapBits(state_sparse.col[index], swap[0], swap[1])
    ordered_state = state_sparse.toarray()
    return ordered_state

#IBMQ.enable_account('267761afd846893dec77bf06dc487d6c7b9569ed20605f88fbddf79e6fb4d149dcfb04da0e4c78994302537adbb9da62c7cfc1d5341f9d4a4c8b416a953bbc9b')
IBMQ.load_account() # Load account from disk
#provider = IBMQ.get_provider(hub='ibm-q-cern')
#backend = provider.get_backend('ibm_cairo')
provider = IBMQ.get_provider(hub='ibm-q')
backend = provider.get_backend('ibm_oslo')
noise_model = NoiseModel.from_backend(backend)
simulator = AerSimulator()
basis_gates = noise_model.basis_gates


N_shots = 1

#Controlled unitary test
"""
qc1 = QuantumCircuit(2)
qc1.x([0,1])
qc1.cx(1, 0)
Z = qc1.to_gate(label = "Z").control(1)

qc = QuantumCircuit(3, 3)
qc.h(0)
qc.append(Z, [0, 1, 2])
qc.measure([0,1,2], [0,1,2])
"""
#Z4 GS Init test
"""
qc = QuantumCircuit(7, 4)
qc.h([2,3])
qc.ccx(1, 3, 0)
qc.cx(2, 0)
qc.cx(3, 1)
qc.measure([0,1,2,3], [0,1,2,3])
"""
#Z3 GS Init test

qc = QuantumCircuit(5, 5)

theta = 2*np.arccos(1/np.sqrt(3))
"""
#fourier right
qc.ry(theta, 2)
qc.ch(2, 3)
qc.x(2)
qc.barrier()

#ground state down
qc.cx(2, 3)
qc.cx(2, 0)
qc.ccx(2, 0, 1)
qc.ccx(2, 1, 0)
qc.cx(3, 0)
qc.ccx(3, 0, 1)
qc.ccx(3, 1, 0)
qc.cx(2, 3)
"""
"""
#Z on both
qc.x(0)
qc.cx(0, 1)
qc.cx(1, 0)
qc.x(2)
qc.cx(2,3)
qc.cx(3,2)
"""

qc.x(0)
qc.x(2)
qc.h(3)
qc.h(4)
qc.z(4)
qc.ccx(0, 2, 4)
qc.ccx(2,1,0)
qc.ccx(3,4,1)


qctrue = qc.copy()

#display(plot_state_hinton(qc))

initial_layout = [1,2,3,4,0]

# Transpile the circuit multiple times for stochastic swap
qc_list = [qc.copy() for i in range(10)]
tqc_list = transpile(qc_list, 
                     simulator, 
                     optimization_level=3, 
                     initial_layout = initial_layout
                     )

# Find best circuit
depths = [tqc.depth() for tqc in tqc_list]
min_depth_index = min(range(len(depths)), key=depths.__getitem__)
tqc = tqc_list[min_depth_index]

display(qc.draw())
tqc.snapshot('final')
display(tqc.draw())

job = simulator.run(tqc, job_name = "simple_test", shots=N_shots)
job_monitor(job)
print("Running time {}s\n".format(job.result().time_taken))
results = job.result()
snap = results.data()['snapshots']['statevector']['final'][0]
#counts = results.get_counts()

#print("\nTotal counts are:", counts)
print("\nTest circuit depth:", qc.depth())
print("Test circuit size:", qc.size())
print("Test circuit width:", qc.width())
print("\nTranspiled circuit depth:", tqc.depth())
print("Transpiled circuit size:", tqc.size())
print("Transpiled circuit width:", tqc.width())


tqctrue = transpile(qctrue, simulator)
tqctrue.snapshot('final')
jobfinal = simulator.run(tqctrue, shots = N_shots)
snaptrue = jobfinal.result().data()['snapshots']['statevector']['final'][0]

print("State backend ", snap)
print("State True", snaptrue)

# Remove small elements from statevector
snap_clean = np.where(abs(snap) > 2.2e-15, snap, 0)

# Convert statevector to sparse matrix
sparse_snap = sp.coo_matrix(snap_clean)

# Unpermute statevector
unpermuted_sparse = ts.swapSparseState(sparse_snap, initial_layout)

# Convert to array
unpermuted_snap = unpermuted_sparse.toarray()[0]

print("Reordered backend ", unpermuted_snap)
print("Fidelity backend and True ", np.abs(np.dot(snap.conjugate(), snaptrue)))
print("Fidelity reordered backend and True ", np.abs(np.dot(unpermuted_snap.conjugate(), snaptrue)))
