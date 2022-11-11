# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 11:12:30 2022

@author: rikci
"""

import scipy.sparse as sp
import tools as ts

N_states_ensemble = 1
ensemble = [sp.load_npz('C:/Users/rikci/.spyder-py3/TESI/SciPy States/test/'+str(i+1)+'.npz')
            for i in range(N_states_ensemble)]
state_ideal = sp.load_npz('C:/Users/rikci/.spyder-py3/TESI/SciPy States/test/ideal.npz')

fid, err = ts.fidelityError(ensemble, state_ideal, 1)
print("\nFidelity with exact state is:", fid, "pm", err)