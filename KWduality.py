# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 18:31:08 2022

@author: rikci
"""

import ToricCode as tc
import tools as ts
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, transpile
from qiskit_ibm_provider import IBMProvider
from qiskit.providers.aer import AerSimulator
from qiskit.tools import job_monitor
from IPython.display import display

provider = IBMProvider()
backend_to_transpile = provider.get_backend('ibmq_manila')
simulator = AerSimulator()

N = 5
qc = QuantumCircuit(2*N, 2*N)
#Initialize odd qubits in + state
[qc.h(2*i+1) for i in range(N)]

#Prepare arbitrary state on even qubits: choose + for simplicity
[qc.h(2*i) for i in range(N)]

# Apply U_CZ cluster state entangler
[qc.cz(2*i, 2*i+1) for i in range(N)]
[qc.cz(2*i+1, 2*i+2) for i in range(N-1)]
qc.cz(2*N-1, 0) # PBC

# Measure 2j on the X basis
[qc.h(2*i) for i in range(N)]
[qc.measure(2*i, 2*i) for i in range(N)]

# Apply string of X on 2j+1 to correct errors
for i in range(N):
    with qc.if_test((2*i, 1)):
        for j in range(i, N-1):
            qc.x(2*j+1)
#     qc.x(1)
#     qc.x(3)
# with qc.if_test((2, 1)):
#     qc.x(3)
    
# Reset 2j qubits
[qc.reset(2*i) for i in range(N)]

qc.snapshot('final')
qc.measure([i for i in range(2*N)], [i for i in range(2*N)])
    

# reg_index = 0
# for x in range(test.nlinks_x):
#     for y in range(test.nlinks_y):
#         control_ancilla = reg_list[reg_index]
#         test.circuit.h(control_ancilla)
#         test.Bp((x,y), power=1, control_qubit=control_ancilla)
#         test.circuit.h(control_ancilla)
#         reg_index += 1

# test.circuit.measure(reg_list, [i for i in range(N_reg)])
display(qc.draw('mpl'))


# Transpile the circuit multiple times for stochastic swap
qc_list = [qc.copy() for i in range(10)]
tqc_list = transpile(qc_list, 
                     simulator, 
                     #optimization_level=2, 
                     )

# Find best circuit, by first removing the circuits that use faulty qubits and
# then choosing the one with smallest depth
"""
# Check faulty cnots (the ones with error probability 1)
props = backend.properties()
cnots = [gate for gate in props.to_dict()['gates'] if gate['gate'] == 'cx']
faulty_cnots = [cnot for cnot in cnots if cnot['parameters'][0]['value'] == 1]
faulty_cnots_qubits = [faulty_cnots[i]['qubits'] for i in range(len(faulty_cnots))]

# Remove circuits that use faulty cnots
clean_tqc_list = []
for tqc in tqc_list:
    faulty = False
    cnots = [instruction for instruction in tqc.data 
             if instruction[0].name == 'cx']
    for op, qubits, clbits in cnots:
        pair = [qubits[0].index, qubits[1].index]
        if pair in faulty_cnots_qubits:
            faulty = True
            break
    if faulty == False: clean_tqc_list.append(tqc)
tqc_list = clean_tqc_list
"""

# Choose circuit of min depth
depths = [tqc.depth() for tqc in tqc_list]
min_depth_index = min(range(len(depths)), key=depths.__getitem__)
tqc = tqc_list[min_depth_index]
display(tqc.draw('mpl'))


# Print circuit depth, size, width before and after transpile
print("\nTest circuit depth:", qc.depth())
print("Test circuit size:", qc.size())
print("Test circuit width:", qc.width())
print("\nTranspiled circuit depth:", tqc.depth())
print("Transpiled circuit size:", tqc.size())
print("Transpiled circuit width:", tqc.width())

# Run circuit on simulator/backend and get results
N_shots = 1000
job = simulator.run(tqc, job_name = "KW Test", shots = N_shots, 
                    #noise_model = noise_model
                    )
job_monitor(job)
result = job.result()
first_ten_states = [result.data()['snapshots']['statevector']['final'][i] for i in range(10)]


# Extract counts
counts = result.get_counts()
print("\nTotal counts are:", counts)
ts.plot_counts(counts)

# Print time taken
print("\nRunning time {}s".format(result.time_taken))



