#hello

import galsim
import numpy as np
import matplotlib.pyplot as plt
import desc.slrealizer
import math

def make_catalog(currLens, currObs):
    image = plot_all_objects(currLens, currObs)

def plot_all_objects(currLens, currObs, save_dir=None):
    MJD, filter, PSF_HWHM, sky_mag = currObs[0], currObs[1], currObs[2], currObs[3]
    filter_AB_offset = 0.00
    curr_galaxy_mag = currLens[currObs[1]+'_SDSS_lens'][0]
    #print 'galaxy magnitude: ', curr_galaxy_mag
    galaxy_flux = np.power(2.5, desc.slrealizer.return_zeropoint()-curr_galaxy_mag+filter_AB_offset)
    #Only one of sigma, fwhm, and half_light_radius may be specified for Gaussian, so we just add the values!
    galaxy = galsim.Gaussian(half_light_radius=currLens['REFF'][0], flux=galaxy_flux)
    big_fft_params = galsim.GSParams(maximum_fft_size=10240)
    phi_angle = currLens['PHIE'][0] * galsim.degrees
    galaxy = galaxy.shear(e=currLens['ELLIP'][0], beta=phi_angle)
    for i in xrange(currLens['NIMG']):
        filter_quasar = currLens[currObs[1]+'_SDSS_quasar']
        curr_lens_mag = -2.5*np.log10(abs(currLens['MAG'][0][i]) + filter_AB_offset) + filter_quasar[0] + filter_AB_offset
        mag_ratio = np.power(2.5, desc.slrealizer.return_zeropoint()-curr_lens_mag)
        lens = galsim.Gaussian(flux=mag_ratio, sigma=0.0) # we are going to convolve later, so we can say that the lens is a point source for now
        lens = lens.shift([currLens['XIMG'][0][i],currLens['YIMG'][0][i]])
        galaxy += lens
    psf = galsim.Gaussian(flux=1, sigma=PSF_HWHM)
    galaxy = galsim.Convolve(galaxy, psf, gsparams=big_fft_params)
    img = galaxy.drawImage(scale=0.2)
    plt.imshow(img.array, interpolation='none', extent=[-10, 10, -10, 10])
    if save_dir is not None:
        plt.savefig(save_dir+'before_deblend.png')
    return img

def generate_data(currLens, currObs, manual_error=True):
    img = desc.slrealizer.plot_all_objects(currLens, currObs)
    MJD, filter, PSF_HWHM, sky_mag = currObs[0], currObs[1], currObs[2], currObs[3]
    RA, RA_err, DEC, DEC_err, first_moment_x_err_calc, first_moment_y_err_calc, size_err = 0, 0, 0, 0, 0, 0, 0
    flux_err_calc = math.pow(10, (22.5 - sky_mag)/2.5)/5 # because Fb = 5 \sigma_b
    try:
        shape_info = img.FindAdaptiveMom()
    except:
        return
    first_moment_x, first_moment_y = (shape_info.moments_centroid.x-(len(img.array)/2.0))*0.2, (shape_info.moments_centroid.y-(len(img.array)/2.0))*0.2 # calculating the real position from the arbitrary pixel position
    # total image intensity for best-fit elliptical Gaussian from adaptive moments.
    flux = shape_info.moments_amp
    size = shape_info.moments_sigma * 0.2 # unit of pixels is returned, so change units for arcsec
    lensID = currLens['LENSID'][0]
    e1 = shape_info.observed_shape.e1
    e2 = shape_info.observed_shape.e2
    e_squared = (e1*e1+e2*e2)
    e = np.sqrt(e_squared)
    if e1 is 0:
        e1 = 0.001 # to solve the division error
    phi = np.arctan(e2/e1)/2.0
    size_err = 0
    if manual_error:
        size += size * noissify_data(desc.slrealizer.get_second_moment_err(), desc.slrealizer.get_second_moment_err_std())
        first_moment_x += noissify_data(desc.slrealizer.get_first_moment_err(), desc.slrealizer.get_first_moment_err_std()) * first_moment_x
        first_moment_y += noissify_data(desc.slrealizer.get_first_moment_err(), desc.slrealizer.get_first_moment_err_std()) * first_moment_y
        flux += flux*noissify_data(desc.slrealizer.get_flux_err(), desc.slrealizer.get_flux_err_std())
        sample_array =  [MJD, filter, RA, RA_err, DEC, DEC_err, first_moment_x, first_moment_x_err_calc, first_moment_y, first_moment_y_err_calc, flux, flux_err_calc, size, size_err, e1, e2, e, phi, PSF_HWHM, sky_mag, lensID]
    return np.array(sample_array)

def noissify_data(mean, stdev):
    """
    Given a mean and a standard deviation of a measurement, add a noise to the data
    """
    return np.random.normal(loc=mean, scale=stdev)
