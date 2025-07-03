import numpy as np
from scipy.optimize import minimize

# Function to generate Lorentzian convolution
def convG(x, x1, fv, dE):
    return 2*np.sqrt(np.log(2))*np.exp(-4*np.log(2)*np.square((x-x1+dE)/fv))/(fv*np.sqrt(np.pi))

# Function to generate Gaussian convolution
def convL(x, x1, fv, dE):
    return fv/(2*np.pi*(np.square(fv/2)+np.square(x-x1+dE)))

# Function to generate psuedo voigt spectra from stick data
def pvoigt_gen(x, G1, L1, L2, dE, dI, SP, stick_x, stick_y):
    # Define FWHM of Voigt function before and after splitting point
    fv1 = 0.5346*2*L1 + 2*np.sqrt(0.2166*L1**2+G1**2)
    fv2 = 0.5346*2*L2 + 2*np.sqrt(0.2166*L2**2+G1**2)

    # Define eta values before and after splitting point
    Eta1 = 1.36603*(2*L1/fv1)-0.47719*pow(2*L1/fv1,2) + 0.11116*pow(2*L1/fv1,3)
    Eta2 = 1.36603*(2*L2/fv2)-0.47719*pow(2*L2/fv2,2) + 0.11116*pow(2*L2/fv2,3)

    # Define the set of points that are represented by 
    x1 = stick_x
    y1 = stick_y
    gx1 = x1[x1>=SP]
    gx2 = x1[x1<SP]
    gy1 = y1[x1>=SP]
    gy2 = y1[x1<SP]

    x11 = np.repeat(x[:, np.newaxis], gx1.shape[0], axis=1)
    x12 = np.repeat(x[:, np.newaxis], gx2.shape[0], axis=1)
    y11 = gy1*((1-Eta1)*convG(x11, gx1, fv1, dE)+Eta1*convL(x11, gx1, fv1, dE))
    y12 = gy2*((1-Eta2)*convG(x12, gx2, fv2, dE)+Eta2*convL(x12, gx2, fv2, dE))

    y = dI*(np.sum(y11, axis=1) + np.sum(y12, axis=1))
    return y

# Function to fit convolution to stick spectra based on experimental spectra using fixed splitting point
def nls_fit_spc_exp(stick_x, stick_y, exp_x, exp_y, sp):

    # bounds = ([0.5, 0.5, 1, -10, 0.01, 6480], [2, 4, 10, 10, 2.0, 6500])
    # bounds = ((0.5,2), (0.5,4), (1,10), (-10,10), (0.01, 2), (6480,6500))
    bounds = ((0.5, 2), (0.5,4), (1,10), (-10,10), (0.01, 2))
    
    # Function to generate spectrum from stick data
    def func(x):
        # Generate spectra
        y = pvoigt_gen(exp_x, *x, sp, stick_x, stick_y)

        # Return loss
        err = np.linalg.norm(y-exp_y)
        # print(err)
        return err
    

    # err = func(x, 0.5, 1.2, 4.0, -4.5, 0.1, 6489)
    # print(err)
    ini_x = np.array([0.5, 1.46, 2.02, 4.29, 1])
    # mat_x = np.array([0.5, 1.49, 4.84, 3.67, 0.878])
    # res = minimize(func, ini_x, bounds=bounds)
    res = minimize(func, ini_x, bounds=bounds)
    
    
    err = func(res.x)
    # mat_err = func(mat_x)
    # print(f"Final Error: {err}")
    # print(f"Matlab Error: {mat_err}")
    return res