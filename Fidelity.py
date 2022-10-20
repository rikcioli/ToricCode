# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 11:12:30 2022

@author: rikci
"""

import scipy.sparse as sp
import tools as ts

N_states_ensemble = 10
ensemble = [sp.load_npz('C:/Users/rikci/.spyder-py3/TESI/SciPy States/Z4 plaq Cairo noise model/'+str(i+1)+'.npz')
            for i in range(N_states_ensemble)]
state_ideal = sp.load_npz('C:/Users/rikci/.spyder-py3/TESI/SciPy States/Z4 plaq Cairo noise model/ideal.npz')

fid, err = ts.fidelityError(ensemble, state_ideal, 2)
print("\nFidelity with exact state is:", fid, "pm", err)