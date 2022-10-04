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



IBMQ.load_account()
provider = IBMQ.get_provider(hub='ibm-q')
backend = provider.get_backend('ibm_nairobi')
noise_model = noise.NoiseModel.from_backend(backend)

# Initialize toric code and registers
N_reg = 3
test = tc.Z4((2,2), (True,True), N_reg)

GScheck = test.circuit.copy()
test.MagneticGS()
GScopy = test.circuit.copy()
test.circuit.snapshot('final')


"""
ancilla = test.N_qubits-test.N_ancillas
reg_list = [ancilla+i for i in range(N_reg)]
test.circuit.h(reg_list)
"""

#Initialize e and m particles and rotate them
"""
test.Z_string((1,1), (3,1))
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


# Print circuit depth, size, width before and after transpile
print("\nTest circuit depth:", qc_list[0].depth())
print("Test circuit size:", qc_list[0].size())
print("Test circuit width:", qc_list[0].width())
print("\nTranspiled circuit depth:", tqc_list[0].depth())
print("Transpiled circuit size:", tqc_list[0].size())
print("Transpiled circuit width:", tqc_list[0].width())



# Run GS circuit on ideal simulator and save state on disk
N_shots = 1
job = simulator.run(tqc_list, job_name = "Ground State", shots = N_shots)
job_monitor(job)
print("Running time {}s\n".format(job.result().time_taken))
snap1 = job.result().data()['snapshots']
state = snap1['statevector']['final'][0]
#save_state(state, name = 'Z4_3bit_new')


#Run GS circuit copy on noisy simulator and measure average fidelity
N_noisy = 1

GScopy.y(1)
GScopy.snapshot('final')
tGScopy = transpile(GScopy, simulator)
job_noisy = simulator.run(tGScopy, job_name = "Ground State noisy", shots = N_noisy)
job_monitor(job_noisy)
print("Running time {}s\n".format(job_noisy.result().time_taken))
snap2 = job_noisy.result().data()['snapshots']
state_noisy_list = snap2['statevector']['final']
avg_fidelity = 0
#state = load_state(name = 'Z4_3bit_new')

for shot in state_noisy_list:
    fidelity = np.abs(np.dot(state.conjugate(), shot))
    avg_fidelity += fidelity/N_noisy
print("Average fidelity with exact state is:", avg_fidelity)

"""
#Check GS copy by evaluating all plaquette and star eigenvalues 
GScheck.set_statevector(state_noisy_list[0])
ancilla = test.N_qubits-test.N_ancillas
reg_list = [ancilla+i for i in range(N_reg)]
GScheck.reset(reg_list)
GScheck.h(reg_list)
GScheck_list = [GScheck.copy() for i in range(6)]
for n_circuit, GS in enumerate(GScheck_list):
    x = 0
    y = 1
    if n_circuit%3 == 2: x = 1
    if n_circuit%3 == 0: y = 0
    for i in range(N_reg):
        if n_circuit < 3:
            test.Bp((x,y), power=3*(2**i), qc = GS, control_qubit=reg_list[i])
        else:
            test.As((x,y), power=3*(2**i), qc = GS, control_qubit=reg_list[i])
    GS.append(QFT(N_reg, inverse=True), [test.circuit.qubits[reg] for reg in reg_list])
    GS.measure(reg_list, [i for i in range(N_reg)])
transp_GScheck = transpile(GScheck_list, simulator)

N_shots_check = 1000
job_GScheck = simulator.run(transp_GScheck, job_name = "GS Check", shots = N_shots_check)
job_monitor(job_GScheck)
result = job_GScheck.result()
print("\nRunning time {}s".format(result.time_taken))
counts_list = [result.get_counts(i) for i in range(6)]
for count in counts_list:
    plot_counts(count)
"""