# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 11:12:30 2022

@author: rikci
"""

import scipy.sparse as sp
import tools as ts

N_states_ensemble = 1000
ensemble = [sp.load_npz('D:/Fisica/TESI/SciPy States/Z4 2 plaq/6e-4 1e-2 100e3 100e3/'+str(i+1)+'.npz')
            for i in range(N_states_ensemble)]
state_ideal = sp.load_npz('D:/Fisica/TESI/SciPy States/Z4 2 plaq/6e-4 1e-2 100e3 100e3/ideal.npz')

fid, err = ts.fidelityError(ensemble, state_ideal, 10)
print("\nFidelity with exact state is:", fid, "pm", err)