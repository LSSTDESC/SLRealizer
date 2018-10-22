"""
The :mod:`moments` module provides functions related to moment estimation used by the :class:`CosmoDC2Realizer` and :class:`SprinkledRealizer` classes.
"""
from __future__ import absolute_import, division, print_function
import numpy as np
import utils
import gc

def sersic_to_mog(sersic_index, size_major, size_minor, e1, e2, vectorize=True):
    '''
    Interface to the function that turns a Sersic profile into a mixture of Gaussians
    and calculates the second moments
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
        half of the rotation angle, or 0.5*arctan(e1/e2), in rad
    vectorize: Boolean
        whether to vectorize the moment calculation over the Gaussian components and the Sersic profiles
    '''
    if type(size_major) == type(size_minor) == type(e1) == type(e2):
        pass
    elif type(size_major) == type(size_minor) == type(e1) == type(e2):
        if isinstance(size_major, (float, int)):
            pass
        else:
            raise ValueError("Only numeric or array types are allowed.")
    else:
        raise ValueError("Arguments must all be single-valued or arrays.")
    
    if sersic_index==4:
        # MoG parameters from Table 1 of Hogg and Lang 2013
        mog_var = np.array([0.00263, 0.01202, 0.04031, 0.12128, 0.36229, 1.23604])**2.0 # v_m
        mog_weights = np.array([0.01308, 0.12425, 0.63551,  2.22560, 5.63989, 9.81523]) # c_m
    elif sersic_index==1:
        mog_var = np.array([0.12068, 0.32730, 0.68542, 1.28089])**2.0 # v_m
        mog_weights = np.array([0.09733, 1.12804, 4.99846, 5.63632]) # c_m
    else:
        raise ValueError("Only de Vaucouleurs and exponential profiles are supported.")
        
    return calculate_mog_second_moments(size_major=size_major, size_minor=size_minor, e1=e1, e2=e2, mog_var=mog_var, mog_weights=mog_weights, vectorize=vectorize)
    
def calculate_mog_second_moments(size_major, size_minor, e1, e2, mog_var, mog_weights, object_id, vectorize):
    '''
    Calculates the second moments of the mixture of Gaussians that approximates the Sersic
    Parameters
    ==========
    mog_var: NumPy.ndarray
        variances of the Gaussians
    mog_weights: NumPy.ndarray
        weights of the Gaussians
    Note
    ====
    See :function:`sersic_to_mog` for the other keyword arguments.
    Return
    =======
    second_moments: NumPy.ndarray
        an array of second moments Ixx, Ixy, Iyy of the mixture of Gaussians
    '''
    num_gaussians = len(mog_var)
    mog_var = mog_var.reshape(1, -1)
    mog_weights = mog_weights.reshape(1, -1)/np.sum(mog_weights)
    
    # TODO: check isinstance(size_major, np.ndarray) and
    if vectorize:
        e, theta = utils.e1e2_to_ephi(e1=e1, e2=e2)
        # Let n be length of parameter array.
        size_major_sq = size_major.reshape(-1, 1)**2.0 # (n, 1) square of r0_major
        size_minor_sq = size_minor.reshape(-1, 1)**2.0 # (n, 1) square of r0_minor
        theta = theta.reshape(-1, 1) # (n, 1)
        cos_theta = np.cos(theta) # (n, 1)
        sin_theta = np.sin(theta) # (n, 1)
        gc.collect()
        
        Ixx_unrotated = np.dot(size_major_sq, mog_var) # (n, num_gaussians), equals r0_major^2*vm
        Iyy_unrotated = np.dot(size_minor_sq, mog_var) # (n, num_gaussians), equals r0_minor^2*vm
        Ixx = cos_theta**2.0*Ixx_unrotated + sin_theta*Iyy_unrotated # (n, num_gaussians), note broadcasting
        Ixy = cos_theta*sin_theta*(Ixx_unrotated - Iyy_unrotated) # (n, num_gaussians), note broadcasting
        Iyy = sin_theta**2.0*Ixx_unrotated + cos_theta**2.0*Iyy_unrotated # (n, num_gaussians), note broadcasting

        second_moments = np.stack([Ixx, Ixy, Iyy], axis=2) # (n, num_gaussians, 3)
        second_moments = np.sum(second_moments*mog_weights.reshape(1, num_gaussians, 1), axis=1) # (n, 3)
        return second_moments
    else:
        # Explicit for-loop over the Gaussian components
        raise NotImplementedError
        #for weight, var in list(zip(mog_weights, mog_var)):
        