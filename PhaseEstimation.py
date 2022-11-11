# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 18:31:08 2022

@author: rikci
"""

import ToricCodeNewFinal as tc
import tools as ts
import numpy as np
from qiskit import transpile, IBMQ
from qiskit.providers.aer import AerSimulator
import qiskit.providers.aer.noise as noise
from qiskit.circuit.library import QFT
from qiskit.tools import job_monitor
from IPython.display import display


"""
IBMQ.load_account()
provider = IBMQ.get_provider(hub='ibm-q')
backend = provider.get_backend('ibm_nairobi')
"""
#IBMQ.enable_account('267761afd846893dec77bf06dc487d6c7b9569ed20605f88fbddf79e6fb4d149dcfb04da0e4c78994302537adbb9da62c7cfc1d5341f9d4a4c8b416a953bbc9b')
provider = IBMQ.get_provider(hub='ibm-q-cern')
backend = provider.get_backend('ibm_cairo')
noise_model = noise.NoiseModel.from_backend(backend)
simulator = AerSimulator(method = "statevector")

# Initialize toric code and registers
N_reg = 3
test = tc.Z4((3,2), (False, False), N_reg)
test.MagneticGS_2()
ancilla = test.N_qubits-test.N_ancillas
reg_list = [ancilla+i for i in range(N_reg)]
test.circuit.h(reg_list)

#Initialize e and m particles and rotate them
test.Z_string((1,1), (2,1))
test.X_string((1,0), (1,1), power = 3)

for i in range(N_reg):
    test.Bp((0,0), power=3*(2**i), control_qubit=reg_list[i])

#Access non trivial sector and measure it with 't Hooft loop
"""
test.Z_string((0,1), (3,1), power = 2)
test.X_string((1,0), (2,0), control_qubit=ancilla)
test.X_string((1,1), (2,1), control_qubit=ancilla)
"""
#Initialize dyons and exchange them an arbitrary number of times
"""
test.init_dyons(low_charges = (1,1), up_charges = (1,1))
for i in range(N_reg):
    test.exchange_countclock(2**(i), control_qubit = reg_list[i])
"""

test.circuit.append(QFT(N_reg, inverse=True), [test.circuit.qubits[reg] for reg in reg_list])
test.circuit.measure(reg_list, [i for i in range(N_reg)])

display(test.circuit.draw('mpl'))

# Transpile the circuit multiple times for stochastic swap
qc_list = [test.circuit.copy() for i in range(200)]
tqc_list = transpile(qc_list, 
                     backend, 
                     optimization_level=3, 
                     )

# Find best circuit, by first removing the circuits that use faulty qubits and
# then choosing the one with smallest depth


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


# Choose circuit of min depth
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
N_shots = 1000
job = simulator.run(tqc, job_name = "Toric Test", shots = N_shots, 
                    noise_model = ts.build_noise()
                    )
job_monitor(job)
result = job.result()

# Extract counts
counts = result.get_counts()
print("\nTotal counts are:", counts)
ts.plot_counts(counts)

# Print time taken
print("\nRunning time {}s".format(result.time_taken))



