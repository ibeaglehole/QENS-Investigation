import numpy as np

q_0_2 = np.loadtxt('q_0.2.txt', unpack = True)
q_0_6 = np.loadtxt('q_0.6.txt', unpack = True) 
q_1_0 = np.loadtxt('q_1.0.txt', unpack = True) 
q_1_4 = np.loadtxt('q_1.4.txt', unpack = True) 
q_1_8 = np.loadtxt('q_1.8.txt', unpack = True) 
q_2_2 = np.loadtxt('q_2.2.txt', unpack = True) 
q_2_6 = np.loadtxt('q_2.6.txt', unpack = True) 
q_3_0 = np.loadtxt('q_3.0.txt', unpack = True)

def Lorentzian_model(omega, A, gamma, omega_0):
    """
    A model of the relationship between the structure factor, S, and energy transfer, omega.
    
    Args:
        omega (array-like): Energy transfer / meV
        A (float): Dimensionless scaling factor
        gamma (float): Describes the width of the curve / meV
        omega_0 (float): Describes the position of the peak of the curve / meV
        
    Returns:
        (array-like): Structure factor as a function of energy transfer, omega.
    """
    denom = np.pi * gamma * (1 + (((omega - omega_0) / (gamma)) ** 2))
    S = A / denom
    return S

def Chudley_Elliot_model(q, D, l):
    """
    Chudley-Elliot diffusion model is used for modelling the variation in gamma (width of curve) as a function of q, which can be used to determine values for the diffusion coefficient, D, and the random walk step size, l.
    
    Args:
        q (array-like): Wavevectors that QENS measurements have been taken at / Å^-1.
        D (float): Diffusion coefficient / Å^2 meV.
        l (float): Random walk step size / Å.
        
    Return:
        (array-like): Describes the width of the curve in a QENS experiment / meV.
    """
    gamma = ((6 * D) / (l ** 2)) * (1 - ((np.sin(q * l)) / (q * l)))
    return gamma

def random_jump_model(q, D, l):
    """
    Random-Jump diffusion model is a simple model describing the diffusion of a random walk in a material.
    
    Args:
        q (array-like): Wavevectors that QENS measurements have been taken at / Å^-1.
        D (float): Diffusion coefficient / Å^2 meV.
        l (float): Random walk step size / Å.
        
    Return:
        (array-like): Describes the width of the curve in a QENS experiment / meV.
    """
    gamma = (6 * D * (q ** 2)) / (1 + ((l ** 2) * (q ** 2)))
    return gamma

def unit_conv(D):
    """
    Function to convert the units of the diffusion coefficient from Å^2 meV to cm^2 s^-1
    
    Args:
        D (float): Diffusion coefficient in units Å^2 meV
        
    Return:
        (float): Diffusion coefficient in units cm^2 s^-1
    """
    return D * (1.519e-4)









