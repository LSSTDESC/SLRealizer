from slrealize_base import SLRealizer
import utils
import constants
import pandas
import numpy as np
import galsim

class OM10Realizer(SLRealizer):

    def __init__(self, catalog, observation):
        super(OM10Realizer, self).__init__(observation=observation)
        self._catalog = catalog
        self._num_systems = len(self._catalog.sample)
        
    def _get_lens(self, lensID=None, rownum=None):
        if lensID is not None and rownum is not None:
            raise ValueError("Need to define either lensID or rownum, not both.")
        
        if lensID is not None:
            return self._catalog.get_lens(lensID)
        elif rownum is not None:
            return self._catalog.get_lens(self._catalog.sample[i])
    
    def draw_lens_random_date(self, lensID=None, rownum=None, save_dir=None):
        """                                                                                                                   
        Draws a lens with the given lensID under randomly chosen observation conditions

        Keyword arguments:
        - lensID: an integer specifying the ID of the lens system in the OM10 catalog

        Returns:
        - A dictionary of the properties of the observed lens
        # TODO say which properties

        """
        # Randomly select observation ID
        obsID = random.randint(0, self.num_obs)
        # Render an image of the system using GalSim
        img, obj = nd.draw_system(self._get_lens(lensID=lensID, rownum=rownum),\
                                  self._get_observation(obsID),\
                                  save_dir)
        # Draw the lens system at the epoch defined by obsID
        fig, axes = plt.subplots(2, figsize=(5, 10))
        axes[0].imshow(img.array, interpolation='none', aspect='auto')
        axes[0].set_title("TRUE MODEL IMAGE, epochID: %d" %obsID)
        params = self.get_lens_params(self.catalog.get_lens(lensID),\
                                       self.observation.loc[obsID])
        # Define a GalSim object using the properties derived from the aggregate system 
        galaxy = galsim.Gaussian(flux=params['flux'], sigma=params['size'])\
                       .shift(float(params['x']), float(params['y']))\
                       .shear(e1=params['e1'], e2=params['e2'])
        # Draw the emulated image
        img = galaxy.drawImage(scale=utils.get_pixel_scale(), method='no_pixel')
        axes[1].imshow(img.array, interpolation='none', aspect='auto')
        axes[1].set_title("EMULATED IMAGE, epochID %d" %obsID)
        if save_dir is not None:
            plt.savefig(save_dir + 'deblending.png')
        return params
    
    def draw_system(self, lens_info, obs_info, save_dir=None):
        '''
        Draws all objects of the given lens system
        in the given observation conditions using GalSim

        Keyword arguments:
        lens_info -- a row of the OM10 DB
        obs_info -- a row of the observation history df
        save_dir -- directory in which to save the image

        Returns:
        A tuple of
        - a Numpy array of the image
        - a GalSim object of the aggreate system used to render
          the image
        '''
        MJD, band, PSF_FWHM, sky_mag = obs_info['expMJD'], obs_info['filter'], obs_info['FWHMeff'], obs_info['fiveSigmaDepth']
        ###############
        # Lens galaxy #
        ###############
        # Compute flux ratio
        gal_mag = lens_info[band + '_SDSS_lens'][0]
        galaxy_flux = from_mag_to_flux(gal_mag\
                                          - return_zeropoint()\
                                          - get_filter_AB_offset())
        # Define GalSim object, give right ellipticity
        galaxy = galsim.Gaussian(half_light_radius=lens_info['REFF_T'][0], flux=galaxy_flux)\
                       .shear(e=lens_info['ELLIP'][0], beta=lens_info['PHIE'][0] * galsim.degrees)
        #################
        # Lensed quasar #
        #################
        for i in xrange(lens_info['NIMG']):
            # TODO ask how to factor in magnification?
            # Compute flux ratio
            curr_obj_mag = lens_info[band + '_SDSS_quasar'][0]\
                           + from_flux_to_mag(lens_info['MAG'][0][i] + get_filter_AB_offset())\
                           + get_filter_AB_offset()
            obj_flux = from_mag_to_flux(curr_obj_mag - return_zeropoint()) # don't subtract AB offset?
            # Define GalSim object, shift accordingly
            lens = galsim.Gaussian(flux=obj_flux, sigma=0.0)\
                         .shift([lens_info['XIMG'][0][i], lens_info['YIMG'][0][i]])
            galaxy += lens
        big_fft_params = galsim.GSParams(maximum_fft_size=10240)
        psf = galsim.Gaussian(flux=1.0, fwhm=PSF_FWHM)
        galsim_obj = galsim.Convolve(galaxy, psf, gsparams=big_fft_params)
        # TODO make pixel scale globally configurable
        # TODO since precision of HSM depends on pixel resolution, either don't pass in img to HSM OR manually define good resolution
        img = galsim_obj.drawImage(scale=get_pixel_scale())
        if save_dir is not None:
            plt.imshow(img.array, interpolation='none', aspect='auto')
            plt.savefig(save_dir+'before_deblend.png')
            plt.close()
        return img, galsim_obj