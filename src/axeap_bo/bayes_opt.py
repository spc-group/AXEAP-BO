from . import xes_fit, xes_calc
import numpy as np
from skopt import gp_minimize

var_names = ["fdd", "fpd", "gpd", "soc", "sov", "tdq_g", "tdq_e", "dt_g", "dt_e", "ds_g", "ds_e", "m_g", "m_e"]

def simul_opt(el, ox, exp_x, exp_y, varis, bounds, n_calls=50, n_initial_points=1, fdd=100.0, fpd=100.0, gpd=100.0, soc=100.0, sov=100.0, tdq_g=0.0, tdq_e=0.0, 
              dt_g=0.0, dt_e=0.0, ds_g=0.0, ds_e=0.0, m_g=0.0, m_e=0.0, em_start=None, em_end=None, G1=1.0, L1=1.2, L2=None, s_pt=None, norm=False):
    var_idxs = np.zeros((len(varis)))
    for i, var in enumerate(varis):
        var_idxs[i] = var_names.index(var)
    print(var_idxs)

    vals = [fdd, fpd, gpd, soc, sov, tdq_g, tdq_e, dt_g, dt_e, ds_g, ds_e, m_g, m_e]

    def var_bo(x):
        vals_i = vals

        for i, idx in enumerate(var_idxs):
            vals_i[int(idx)] = x[i]
        
        stick_x, stick_y = xes_calc.gen_sticks(el, ox, *vals)
        res = xes_fit.min_fit(stick_x, stick_y, exp_x, exp_y)
        
        spec = xes_fit.pvoigt_gen(exp_x, *res.x, stick_x, stick_y)
        err = np.linalg.norm(spec-exp_y)
        print(vals_i)
        print(err)
        return err
    
    res = gp_minimize(var_bo, bounds, n_calls=n_calls, n_initial_points=n_initial_points)
    print(res)
    return res

def seq_opt(el, ox, exp_x, exp_y, vars_seq, bounds_seq, calls_seq, n_loops=10, fdd=100.0, fpd=100.0, gpd=100.0, soc=100.0, sov=100.0, tdq_g=0.0, tdq_e=0.0, dt_g=0.0, dt_e=0.0, ds_g=0.0, ds_e=0.0,
                m_g=0.0, m_e=0.0, em_start=None, em_end=None, G1=1.0, L1=1.2, L2=None, s_pt=None, norm=False):
    
    vals = [fdd, fpd, gpd, soc, sov, tdq_g, tdq_e, dt_g, dt_e, ds_g, ds_e, m_g, m_e]
    print(calls_seq)
    for j in range(n_loops):
        for i, varis in enumerate(vars_seq):
            # Single optimization step
            bounds = bounds_seq[i]
            calls = calls_seq[i]
            print(calls)
            res = simul_opt(el, ox, exp_x, exp_y, varis, bounds, calls, 1, *vals)
            
            # Store optimization result
            var_idxs = np.zeros((len(varis)))
            for i, var in enumerate(varis):
                var_idxs[i] = var_names.index(var)
            
            for i, idx in enumerate(var_idxs):
                vals[int(idx)] = res.x[i]
            print(vals)

            
