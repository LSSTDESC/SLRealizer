import numpy as np
from fractions import Fraction
import math
import numpy as np
import pandas

def from_flux_to_mag(flux_ratio):
    return (-2.5)**flux_ratio

def from_mag_to_flux(mag_diff):
    return 10.0**(-0.4*mag_diff)

def noissify_data(mean, stdev, measurement=1.):
    """
    Given a mean and a standard deviation of a measurement, add a noise to the data
    """
    return measurement*np.random.normal(loc=mean, scale=stdev)

def get_pixel_scale():
    return 0.2
'''
def return_coordinate(first_moment_x, first_moment_y):
    """
    This method returns the RA and DEC values that take the x and y offset into account.
    RA and DEC values are provided in the observational history,
    and this method also assumes a random error in the measurement (+1/-1 deg)
    """

###                                                                         
    pos_err = 0.0 # unit : degree                                               
    pos_err_std = Fraction(1, 3) # Assuming that the three sigma is one degree, one sigma is 0.3333 degree                                                     
    ###      

    real_coordinate = desc.slrealizer.return_obs_RA_DEC()
    RA = real_coordinate.ra.deg
    DEC = real_coordinate.dec.deg
    # add the offset by the first moment
    RA += first_moment_x/(3600*np.cos(DEC))
    DEC += first_moment_y/3600
    # draw random number to set position in the FOV
    RA += np.random.uniform(-1.75, 1.75)
    DEC += np.random.uniform(-1.75, 1.75)
    RA += noissify_data(pos_err, pos_err_std)
    DEC += noissify_data(pos_err, pos_err_std)
    RA_err = 0
    DEC_err = 0
    return RA, RA_err, DEC, DEC_err
'''
def return_mean_properties(lens_array):
    """
    returns the mean value of flux, 1st moment of x, 1st moment of y, qzz, qxy, qyy, flux_err, 1st moment x error, 1st moment y error, qxx error, qxy error, and qyy error
    """
    # TODO do not hardcode order, use list comprehension for fixed order?
    return lens_array['flux'].mean(), lens_array['x'].mean(), lens_array['y'].mean(), lens_array['size'].mean(), lens_array['flux_err'].mean(), lens_array['x_com_err'].mean(), lens_array['y_com_err'].mean(), lens_array['size_err'].mean(), lens_array['e1'].mean(), lens_array['e2'].mean()
