from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from realize_sl import SLRealizer
from utils.constants import *
from utils.utils import *
import numpy as np
import pandas as pd
import galsim
import gc # need this to optimize memory usage

class OM10Realizer(SLRealizer):
    
    """
    
    A class that realizes objects in the OM10 mock quasar catalog
    under the given observation conditions, 
    into LSST DRP Source and Object catalogs
    
    """
    
    def __init__(self, observation, catalog, debug=False, add_moment_noise=True, add_flux_noise=True):
        #super(OM10Realizer, self).__init__(observation) # Didn't work for some reason
        self.as_super = super(OM10Realizer, self)
        self.as_super.__init__(observation, add_moment_noise=add_moment_noise, add_flux_noise=add_flux_noise)
        self.catalog = catalog
        self.num_systems = len(self.catalog.sample)
        self.DEBUG = debug
        
    def get_lens_info(self, objID=None, rownum=None):
        if objID is not None and rownum is not None:
            raise ValueError("Need to define either objID or rownum, not both.")
        
        if objID is not None:
            return self.catalog['objectId'] # Don't use this (?)
        elif rownum is not None:
            return self.catalog.sample[rownum]
       
    def _om10_to_galsim(self, lens_info, band):
        """
        Converts OM10's column values into GalSim terms
        
        Keyword arguments:
        lens_info -- a row of the OM10 DB
        band -- the filter used to observe

        Returns:
        A dictionary (named galsimInput) containing properties that 
        can be passed into GalSim (See code for which properties)
        """
        
        # We will input flux in units of nMgy
        mag   = lens_info[band + '_SDSS_lens'] 
        flux  = mag_to_flux(mag, to_unit='nMgy') # in nMgy
        hlr   = lens_info['REFF_T'] # REFF_T is in arcsec
        e     = lens_info['ELLIP']
        beta  = lens_info['PHIE'] * galsim.degrees # PHIE is in degrees
        
        galsimInput={'flux': flux,
                       'half_light_radius': hlr,
                       'e': e,
                       'beta': beta,
                       'num_objects': lens_info['NIMG']}
        
        for obj in xrange(lens_info['NIMG']):
            obj_mag = lens_info[band + '_SDSS_quasar']\
                      + flux_to_mag(abs(lens_info['MAG'][obj]))
                      #+ flux_to_mag(lens_info['MAG'][obj] + get_filter_AB_offset())
            galsimInput['flux_'+str(obj)] = mag_to_flux(obj_mag, to_unit='nMgy') # don't subtract AB offset?
            galsimInput['xy_'+str(obj)] = lens_info['XIMG'][obj], lens_info['YIMG'][obj]
            
        return galsimInput
    
    def _om10_to_lsst(self, obs_info, lens_info):
        """Converts OM10 column values into LSST source table format
        using analytical mooment calculation

        Keyword arguments:
        obs_info -- dictionary containing the observation conditions
        lens_info -- dictionary containing the lens properties
        
        Returns:
        A dictionary (named derivedProps) containing properties that 
        can be used to propagate one row of the source table
        """
        
        derivedProps = {}
        
        histID, MJD, band, PSF_FWHM, sky_mag = obs_info
        
        numQuasars = lens_info['NIMG']
        lens_mag = lens_info[band + '_SDSS_lens']
        lens_flux = mag_to_flux(lens_mag, to_unit='nMgy')
        q_mag_arr = lens_info[band + '_SDSS_quasar'] + flux_to_mag(np.abs(np.array(lens_info['MAG'][:numQuasars])))
        q_flux_arr = mag_to_flux(q_mag_arr, to_unit='nMgy')
        q_tot_flux = np.sum(q_flux_arr)
        derivedProps['apFlux'] = lens_flux + q_tot_flux
        derivedProps['apMag'] = flux_to_mag(derivedProps['apFlux'], from_unit='nMgy')
        
        #TODO
        derivedProps['apMag'] = 0.0
        derivedProps['apMagErr'] = 0.0
        
        #################################
        # Analytical moment calculation #
        #################################
        
        # Flux ratios weight the moment contributions from the lens and each quasar image
        lensFluxRatio = lens_flux/derivedProps['apFlux']
        qFluxRatios = q_flux_arr/derivedProps['apFlux']
        
        # Slice x, y positions to number of quasars
        lens_info['XIMG'] = lens_info['XIMG'][:numQuasars]
        lens_info['YIMG'] = lens_info['YIMG'][:numQuasars]
        
        # Convert Gaussian FWHM, HLR to sigma of the covariance matrix
        sigmasq_psf = fwhm_to_sigma(PSF_FWHM)**2.0
        sigmaSqLens = hlr_to_sigma(lens_info['REFF_T'])**2.0
        
        # Compute first moment contributions from the lens and each quasar image
        lensIxContrib = lensIyContrib = 0.0 # because lens is centered at (0, 0)
        qIxContrib = qFluxRatios*lens_info['XIMG']
        qIyContrib = qFluxRatios*lens_info['YIMG']
        derivedProps['x'] = lensIxContrib + np.sum(qIxContrib)
        derivedProps['y'] = lensIyContrib + np.sum(qIyContrib)
        
        # Compute second moment contributions from the lens and each quasar image
        lensIxxContrib = lensFluxRatio*(sigmasq_psf + sigmaSqLens + derivedProps['x']**2.0)
        lensIyyContrib = lensFluxRatio*(sigmasq_psf + sigmaSqLens + derivedProps['y']**2.0)
        lensIxyContrib = lensFluxRatio*derivedProps['x']*derivedProps['y'] # because lens is centered at (0, 0)
        qIxxContrib = qFluxRatios*((lens_info['XIMG'] - derivedProps['x'])**2.0 + sigmasq_psf)
        qIyyContrib = qFluxRatios*((lens_info['YIMG'] - derivedProps['y'])**2.0 + sigmasq_psf)
        qIxyContrib = qFluxRatios*(lens_info['XIMG'] - derivedProps['x'])*(lens_info['YIMG'] - derivedProps['y'])
        
        # Add to get total second moments
        Ixx = lensIxxContrib + np.sum(qIxxContrib)
        Iyy = lensIyyContrib + np.sum(qIyyContrib)
        Ixy = lensIxyContrib + np.sum(qIxyContrib)
        derivedProps['trace'] = Ixx + Iyy
        if self.DEBUG: derivedProps['det'] = Ixx*Iyy - Ixy**2.0
        
        # Compute ellipticities from second moments
        denom = Ixx + Iyy #+ 2.0*(Ixx*Iyy - Ixy**2.0)**0.5 # uncommenting gives g1, g2 
        derivedProps['e1'] = (Ixx - Iyy)/denom
        derivedProps['e2'] = 2.0*Ixy/denom
        
        if self.DEBUG:
            derivedProps['Ixx'] = Ixx
            derivedProps['Iyy'] = Iyy
            derivedProps['Ixy'] = Ixy
        
        return derivedProps
            
    def draw_system(self, obs_info, lens_info, save_dir=None):
        galsimInput = self._om10_to_galsim(lens_info, obs_info['filter'])
        return self.as_super.draw_system(galsimInput=galsimInput, obs_info=obs_info, save_dir=save_dir)
    
    def estimate_hsm(self, obs_info, lens_info):
        """
        Performs GalSim's HSM shape estimation on the image
        rendered with lens properties in lens_info
        under the observation conditions in obs_info
        
        Keyword arguments:
        obs_info -- dictionary containing the observation conditions
        lens_info -- dictionary containing the lens properties 
        
        Returns
        a dictionary containing the shape information 
        numerically derived by HSM
        """
        galsim_img = self.draw_system(lens_info=lens_info, obs_info=obs_info, save_dir=None)
        return self.as_super.estimate_hsm(galsim_img=galsim_img)
    
    def draw_emulated_system(self, obs_info, lens_info):
        """
        Draws the emulated system, i.e. draws the aggregate system
        from properties HSM derived from the image, which was in turn
        drawn from the catalog's truth properties.
        Only runs when DEBUG == True.
        
        Returns
        a GalSim Image object of the emulated system
        """
        hsmOutput = self.estimate_hsm(obs_info, lens_info)
        return self.as_super.draw_emulated_system(hsmOutput)
    
    def create_source_row(self, obs_info, lens_info, use_hsm=False):
        '''
        Returns a dictionary of lens system's properties
        computed the image of one lens system and the observation conditions,
        which makes up a row of the source table.

        Keyword arguments:
        image -- a Numpy array of the lens system's image
        obs_info -- a row of the observation history df
        
        Returns
        A dictionary with properties derived from HSM estimation
        (See code for which properties)
        '''
        objectId = lens_info['LENSID']
        if use_hsm:
            derivedProps = self.estimate_hsm(obs_info, lens_info)
            if derivedProps is None:
                return None
        else:
            derivedProps = self._om10_to_lsst(obs_info=obs_info, lens_info=lens_info)
        
        return self.as_super.create_source_row(derivedProps=derivedProps, objectId=objectId, obs_info=obs_info)
    
    def make_source_table_vectorized(self, output_source_path, include_time_variability):
        """
        Generates the source table and saves it as a csv file.

        Keyword arguments:
        output_source_path -- save path for the output source table
        include_time_variability -- whether to include intrinsic quasar variability
        
        Returns (only if self.DEBUG == True):
        a Pandas dataframe of the source table
        """
        import time
        
        start = time.time()
        self._preformat_source_table()
        
        if include_time_variability:
            self.include_quasar_variability(save_output=False)
        self.source_table.reset_index(inplace=True)
        src = self.source_table
        
        # Convert each quasar magnitude back to flux
        for q in range(4):
            src['q_flux_' + str(q)] = mag_to_flux(src['q_mag_' + str(q)], to_unit='nMgy')
        # Get total flux
        q_flux_cols = ['q_flux_' + str(q) for q in range(4)]
        src['apFlux'] = src[q_flux_cols + ['lens_flux']].sum(axis=1)
        
        self.source_table = src
        self._include_moments()
        
        # Add flux noise
        src['apFluxErr'] = mag_to_flux(src['fiveSigmaDepth']-22.5)/5.0 # because Fb = 5 \sigma_b
        if self.add_flux_noise:
            src['apFlux'] += add_noise(mean=0.0, 
                                       stdev=src['apFluxErr'], 
                                       shape=src['apFluxErr'].shape)
        # Get total magnitude
        src['apMag'] = flux_to_mag(src['apFlux'], from_unit='nMgy')
        # Propagate to get error on magnitude
        src['apMagErr'] = (2.5/np.log(10.0)) * src['apFluxErr'] / src['apFlux']
        gc.collect()
            
        # Remove remaining unused columns
        src.drop(['lensFluxRatio', 'lens_mag', 'q_mag', 'sigmasq_lens', 'sigmasq_psf', 'NIMG',], axis=1, inplace=True)
        for q in range(4):
            src.drop(['MAG_%d'%q] + ['XIMG_%d'%q] + ['YIMG_%d'%q] + ['qFluxRatio_%d'%q], axis=1, inplace=True)
        gc.collect()

        ############################################
        # Final column reordering & saving to file #
        ############################################
        src = src[self.source_columns]
        if self.DEBUG:
            out_num_lenses = src['objectId'].nunique()
            out_num_times = src['MJD'].nunique()
            print("Result of making source table: ")
            print("Number of observations: ", out_num_times)
            print("Number of lenses: ", out_num_lenses)
        src.set_index('objectId', inplace=True)
        
        src.to_csv(output_source_path)
        gc.collect()
        end = time.time()

        print("Done making the source table with %d row(s) in %0.2f seconds using vectorization." %(len(src), end-start))
        self.sourceTable = src
        if self.DEBUG:
            return src
    
    def _preformat_source_table(self):
        """
        Initializes self.source_table with the column conventions
        that can be used by SLRealizer's helper functions
        """
        from astropy.table import Table, Column
        
        catalogAstropy = self.catalog.sample # the Astropy table underlying OM10 object
        
        #################################
        # OM10 --> Pandas DF conversion #
        #################################
        lensMagCols = [b + '_SDSS_lens' for b in 'ugriz']
        qMagCols = [b + '_SDSS_quasar' for b in 'ugriz']
        saveCols = lensMagCols + qMagCols + ['REFF_T', 'NIMG', 'LENSID', 'ELLIP', 'PHIE']
        saveValues = [catalogAstropy[c] for c in saveCols]
        saveColDict = dict(zip(saveCols, saveValues))
        collapsedColDict = get_1D_columns(multidimColNames=['MAG', 'XIMG', 'YIMG'], table=catalogAstropy)
        saveColDict.update(collapsedColDict)
        catalog = Table(saveColDict.values(), names=saveColDict.keys()).to_pandas()
        catalog.drop_duplicates('LENSID', inplace=True)

        ####################################
        # Merging catalog with observation #
        ####################################
        observation = self.observation.copy()
        catalog['key'] = 0
        observation['key'] = 0
        src = catalog.merge(observation, how='left', on='key')
        src.drop('key', 1, inplace=True)
        gc.collect()
        
        ##############################################
        # Rename columns and set null values to zero #
        ##############################################
        #src.drop(['fiveSigmaDepth', ], axis=1, inplace=True)
        src.rename(columns={'obsHistID': 'ccdVisitId',
                            'LENSID': 'objectId',
                            'expMJD': 'MJD',
                            'FWHMeff': 'psf_fwhm',
                            'ELLIP': 'e',
                            'PHIE': 'beta',
                             }, inplace=True)
        gc.collect()
        
        # Set unused band magnitudes to zero, to
        # work only with the magnitude in the observed filter
        for b in 'ugriz': # b = observed filter
            setZeroLens = lensMagCols[:]
            setZeroLens.remove(b + '_SDSS_lens')
            setZeroQ = qMagCols[:]
            setZeroQ.remove(b + '_SDSS_quasar')
            src.loc[src['filter'] == b, setZeroLens] = 0.0
            src.loc[src['filter'] == b, setZeroQ] = 0.0
        src['lens_mag'] = src[lensMagCols].sum(axis=1)
        src['q_mag'] = src[qMagCols].sum(axis=1)
        src.drop(lensMagCols + qMagCols, axis=1, inplace=True)
        gc.collect()
        
        # Convert magnitudes into fluxes
        src['lens_flux'] = mag_to_flux(src['lens_mag'], to_unit='nMgy')
        for q in range(4):
            src['q_mag_' + str(q)] = src['q_mag'] + flux_to_mag(np.abs(src['MAG_' + str(q)]))
            src['q_flux_' + str(q)] = mag_to_flux(src['q_mag_' + str(q)], to_unit='nMgy')

        # Set fluxes of nonexistent quasar images to zero
        src.loc[src['NIMG'] == 2, ['q_flux_2', 'q_flux_3']] = 0.0
        src.loc[src['NIMG'] == 3, ['q_flux_3']] = 0.0
        
        self.source_table = src
    
        
    
    #def add_time_variability INHERITED
    #def make_source_table INHERITED
    #def compare_truth_vs_emulated INHERITED
