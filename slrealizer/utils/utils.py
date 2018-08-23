import numpy as np


def get_first_moments_from_image(image_array, pixel_scale):
    """ 
    Returns the first moments in arcsec units numerically computed from
    an image on a pixel grid
    
    Keyword arguments:
    image_array -- a numpy image array
    pixel_scale -- scale factor for the image in arcsec/pixel
 
    Returns:
    a tuple of the first moments Ix, Iy in arcsec
    """
    nx, ny = image_array.shape
    offset = np.array([(nx - 1)/2, (ny - 1)/2]).reshape(2, 1, 1)
    pixel_grid = np.indices([nx, ny]) - offset # shape (2, nx, ny)
    x_coords = pixel_grid[1, :, :] * pixel_scale
    y_coords = pixel_grid[0, :, :] * pixel_scale
    total_flux = np.sum(image_array)
    Ix = np.sum(image_array * x_coords) / total_flux 
    Iy = np.sum(image_array * y_coords) / total_flux 
    return Ix, Iy

def get_second_moments_from_image(image_array, pixel_scale):
    """ 
    Returns the second moments in arcsec units numerically computed from
    an image on a pixel grid
    
    Keyword arguments:
    image_array -- a numpy image array
    pixel_scale -- scale factor for the image in arcsec/pixel
 
    Returns:
    a tuple of the second moments Ixx, Ixy, Iyy in arcsec
    """
    nx, ny = image_array.shape
    Ix, Iy = get_first_moments_from_image(image_array, pixel_scale=1.0)
    offset = np.array([(ny - 1)/2 + Iy, (nx - 1)/2 + Ix]).reshape(2, 1, 1)
    pixel_grid = np.indices([ny, nx]) - offset # shape (2, nx, ny)
    x_coords = pixel_grid[1, :, :] * pixel_scale
    y_coords = pixel_grid[0, :, :] * pixel_scale
    total_flux = np.sum(image_array)
    Ixx = np.sum(image_array * np.power(x_coords, 2.0)) / total_flux 
    Ixy = np.sum(image_array * x_coords * y_coords) / total_flux 
    Iyy = np.sum(image_array * np.power(y_coords, 2.0)) / total_flux 
    return Ixx, Ixy, Iyy

def e1e2_to_ephi(e1, e2):
    e = np.power(np.power(e1, 2.0) + np.power(e2, 2.0), 0.5)
    phi = 0.5*np.arctan(e2/e1)
    return e, phi

def ephi_to_e1e2(e, phi):
    e1 = e*np.cos(2.0*phi)
    e2 = e*np.sin(2.0*phi)
    return e1, e2

def get_1D_columns(multidimColNames, table):
    totalColDict = {}
    for mc in multidimColNames:
        colSize = 4 #table[mc].shape[1] is always 4
        colNames = [mc + '_' + str(c) for c in range(colSize)]
        colValues = [table[mc][:, c].data for c in range(colSize)]
        colDict = dict(zip(colNames, colValues))
        totalColDict.update(colDict)
    return totalColDict

def hlr_to_sigma(hlr):
    return hlr/np.sqrt(2.0*np.log(2.0))

def fwhm_to_sigma(fwhm):
    return fwhm/np.sqrt(8.0*np.log(2.0))

def pixel_to_physical(pixelPos, canvas_size, pixel_scale):
    return (pixelPos - 0.5*canvas_size - 0.5)*pixel_scale

def physical_to_pixel(physicalPos, canvas_size, pixel_scale):
    return physicalPos/pixel_scale + 0.5*canvas_size + 0.5

def scale_mag_as_flux(mag, flux_scale=1.0):
    """
    Identical to from_flux_to_mag(from_mag_to_flux(mag)*flux_scale)
    """
    return mag - 2.5*np.log10(flux_scale)

def flux_to_mag(flux, zeropoint_mag=0.0, from_unit=None, to_unit=None):
    if from_unit=='nMgy':
        zeropoint_mag=22.5
    return zeropoint_mag-2.5*np.log10(flux)

def mag_to_flux(mag, zeropoint_mag=0.0, from_unit=None, to_unit=None):
    if to_unit=='nMgy':
        zeropoint_mag=22.5
    return np.power(10.0, -0.4*(mag - zeropoint_mag))

def add_noise(mean, stdev, shape=None, measurement=1.0):
    """
    Given a mean and a standard deviation of a measurement, adds Gaussian noise to the data
    
    Keyword arguments:
    mean -- the mean of Gaussian
    stdev -- the standard deviation of Gaussian
    shape -- the array shape of noise to be returned
             If None, returns a scalar noise [default: None]
    measurement -- scaling factor, for adding fractional errors. 
                   If 1.0, error is absolute. [default: 1.0]
    """
    return measurement*np.random.normal(loc=mean, scale=stdev, size=shape)

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
