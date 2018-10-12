"""
The :mod:`moments` module provides functions related to moment estimation used by the :class:`CosmoDC2Realizer` and :class:`SprinkledRealizer` classes.
"""
from __future__ import absolute_import, division, print_function
import numpy as np
import gc

def sersic_to_mog(sersic_index=4, size_major, size_minor, e, theta):
    '''
    Interface to functions that turn a Sersic profile into a mixture of Gaussians
    Parameters
    ==========
    sersic_index: int or NumPy.ndarray of int
        the Sersic index indicating type of Sersic profile (1 for exponential, 4 for de Vaucouleurs)
    size_major: float or NumPy.ndarray of float
        the half-light radius of major axis
    size_minor: float or NumPy.ndarray of float
        the half-light radius of minor axis
    e: float or NumPy.ndarray of float
        ellipticity modulus
    theta: float or NumPy.ndarray of float
        half of the rotation angle, or 0.5*arctan(e1/e2)
    Returns
    =======
    second_moments: NumPy.ndarray
        an array of second moments Ixx, Ixy, Iyy of the mixture of Gaussians
    '''
    if type(size_major) == type(size_minor) == type(e) == type(theta) == np.ndarray:
        pass
    elif type(size_major) == type(size_minor) == type(e) == type(theta):
        if isinstance(size_major, (float, int)):
            pass
        else:
            raise ValueError("Only numeric or array types are allowed.")
    else:
        raise ValueError("Arguments must all be single-valued or arrays.")
    
    if sersic_index==4:
        return devaucouleurs_to_mog(size_major=size_major, size_minor=size_minor, e1=e1, e2=e2)
    elif sersic_index==1:
        return exponential_to_mog(size_major=size_major, size_minor=size_minor, e1=e1, e2=e2)
    else:
        raise ValueError("Only de Vaucouleurs and exponential profiles are supported.")
    
def devaucouleurs_to_mog(size_major, size_minor, e1, e2, vectorize=True):
    
    # MoG parameters from Table 1 of Hogg and Lang 2013
    dev_mog_stds = np.array([0.00263, 0.01202, 0.04031, 0.12128, 0.36229, 1.23604]).reshape(1, -1) # sqrt(v_m)
    dev_mog_weights = np.array([0.01308, 0.12425, 0.63551,  2.22560, 5.63989, 9.81523]).reshape(1, -1) # c_m
    num_gaussians = 6
    
    if isinstance(size_major, np.ndarray):
        if vectorize:
            size_major = size_major.reshape(-1, 1)**2.0 # (n, 1) if n is length of parameter array
            gaussian_second_moments = np.dot(size_major, dev_mog_stds**2.0) # (n, 6)
            
            
        
        
        