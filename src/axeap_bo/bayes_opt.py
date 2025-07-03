from . import xes_fit, xes_calc
import numpy as np
from skopt import gp_minimize

var_names = ["fdd", "fpd", "gpd", "soc", "sov", "tdq_g", "tdq_e", "dt_g", "dt_e", "ds_g", "ds_e", "m_g", "m_e"]

def simul_opt(vars, bounds, el, ox, fdd=100.0, fpd=100.0, gpd=100.0, soc=100.0, sov=100.0, tdq_g=0.0, tdq_e=0.0, dt_g=0.0, dt_e=0.0, ds_g=0.0, ds_e=0.0,
                m_g=0.0, m_e=0.0, em_start=None, em_end=None, G1=1.0, L1=1.2, L2=None, s_pt=None, norm=False):
    var_idxs = np.zeros((len(vars)))
    for i, var in enumerate(vars):
        var_idxs[i] = var_names.index(var)
    print(var_idxs)

    vals = [fdd, fpd, gpd, soc, sov, tdq_g, tdq_e, dt_g, dt_e, ds_g, ds_e, m_g, m_e]

    def var_bo(x):
        for i, idx in enumerate(var_idxs):
            print(i, idx)
            print(vals)
            print(x)
            vals[int(idx)] = x[i]
        
        stick_x, stick_y, spec_x, spec_y, = simul_opt(el, ox, *vals)
        return 0.5
    
    res = gp_minimize(var_bo, bounds, n_calls=10)
    return res

def seq_opt(vars_seq, bounds_seq, el, ox, fdd=100.0, fpd=100.0, gpd=100.0, soc=100.0, sov=100.0, tdq_g=0.0, tdq_e=0.0, dt_g=0.0, dt_e=0.0, ds_g=0.0, ds_e=0.0,
                m_g=0.0, m_e=0.0, em_start=None, em_end=None, G1=1.0, L1=1.2, L2=None, s_pt=None, norm=False):
    
    for i, vars in enumerate(vars_seq):
        bounds = bounds_seq[i]
        res = simul_opt(vars, bounds, el, ox, fdd=100.0, fpd=100.0, gpd=100.0, soc=100.0, sov=100.0, tdq_g=0.0, tdq_e=0.0, dt_g=0.0, dt_e=0.0, ds_g=0.0, ds_e=0.0,
                m_g=0.0, m_e=0.0, em_start=None, em_end=None, G1=1.0, L1=1.2, L2=None, s_pt=None, norm=False)
        
        print(res)

            
