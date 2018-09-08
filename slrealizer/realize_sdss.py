from __future__ import absolute_import, division, print_function

from slrealizer import SLRealizer
import slrealizer.utils.constants as constants
import slrealizer.utils.utils as utils
import pandas as pd
import numpy as np

class SDSSRealizer(SLRealizer):

    """
    A class that realizes objects in the SDSS DR catalog
    under the given observation conditions,
    into LSST DRP Source and Object catalogs.
    """

    def __init__(self, observation, catalog, debug=False, add_moment_noise=True, add_flux_noise=True):
        #super(SDSSRealizer, self).__init__(observation) # Didn't work for some reason
        self.as_super = super(SDSSRealizer, self)
        self.as_super.__init__(observation, add_moment_noise=add_moment_noise, add_flux_noise=add_flux_noise)
        self.catalog = catalog
        self.num_systems = len(self.catalog)
        self.DEBUG = debug
        self.sdss_pixel_scale = 0.396 #arcsec/pixel

    def get_lens_info(self, objID=None, rownum=None):
        if objID is not None and rownum is not None:
            raise ValueError("Need to define either objID or rownum, not both.")

        if objID is not None:
            return self.catalog.get_lens(objID) # Don't use this (?)
        elif rownum is not None:
            return self.catalog.loc[rownum]

    def _sdss_to_galsim(self, lens_info, band):
        raise NotImplementedError

    def _sdss_to_lsst(self, obs_info, lens_info):
        raise NotImplementedError

    def draw_system(self, obs_info, lens_info, save_dir=None):
        raise NotImplementedError

    def estimate_hsm(self, obs_info, lens_info):
        raise NotImplementedError

    def draw_emulated_system(self, obs_info, lens_info):
        raise NotImplementedError

    def create_source_row(self, obs_info, lens_info, method="analytical"):

        histID, MJD, band, PSF_FWHM, sky_mag = obs_info
        row = {}

        row['filter'] = band
        row['ccdVisitId'] = histID
        row['MJD'] = MJD
        row['psf_fwhm'] = PSF_FWHM
        row['apFluxErr'] = utils.mag_to_flux(sky_mag-22.5) / 5.0
        row['objectId'] = lens_info['objectId']
        row['apFlux'] = lens_info['modelFlux_' + band]
        row['x'] = np.cos(np.deg2rad(lens_info['offsetDec_' + band]*3600.0)) * lens_info['offsetRa_' + band]
        row['y'] = lens_info['offsetDec_' + band]
        row['trace'] = lens_info['mRrCc_' + band]*(self.sdss_pixel_scale)**2.0 + 2.0*utils.fwhm_to_sigma(PSF_FWHM)**2.0
        row['e1'] = lens_info['mE1_' + band]
        row['e2'] = lens_info['mE2_' + band]

        return row

    def make_source_table_vectorized(self, save_file):
        import gc # need this to optimize memory usage
        import time

        start = time.time()

        ####################################
        # Merging catalog with observation #
        ####################################
        catalog = self.catalog.copy()
        observation = self.observation.copy()
        catalog['key'] = 0
        observation['key'] = 0
        src = catalog.merge(observation, how='left', on='key')
        src.drop('key', 1, inplace=True)
        gc.collect()

        ####################################
        # Collapsing multi-band properties #
        # into one of observed band        #
        ####################################
        setZeroDict = {}
        propsToCollapse = ['modelFlux', 'offsetRa', 'offsetDec', 'mRrCc', 'mE1', 'mE2', ]
        # Initialize dictionary of columns we want to collapse
        for p in propsToCollapse:
            setZeroDict[p] = [p + '_' + b for b in 'ugriz']
        # Set unused column values to zero
        for b in 'ugriz': # b = observed filter
            for p in propsToCollapse: # multi-band columns to collapse
                setZeroCols = setZeroDict[p][:]
                setZeroCols.remove(p + '_' + b)
                src.loc[src['filter'] == b, setZeroCols] = 0.0
        # Collapse
        for p in propsToCollapse: # multi-band columns to collapse
            src[p] = src[[p + '_' + b for b in 'ugriz']].sum(axis=1)
            src.drop([p + '_' + b for b in 'ugriz'], axis=1, inplace=True)
        gc.collect()

        ################
        # Adding noise #
        ################
        src['apFluxErr'] = utils.mag_to_flux(src['fiveSigmaDepth'] - 22.5)/5.0
        if self.add_flux_noise:
            src['modelFlux'] += utils.add_noise(mean=0.0,
                                          stdev=src['apFluxErr'],
                                          shape=src['apFluxErr'].shape) # flux rms not skyEr
        src['x'] = np.cos(np.deg2rad(src['offsetDec']*3600.0))*src['offsetRa']
        src['y'] = src['offsetDec']
        src['trace'] = src['mRrCc']*(self.sdss_pixel_scale**2.0) + 2.0*np.power(utils.fwhm_to_sigma(src['FWHMeff']), 2.0)
        if self.add_moment_noise:
            src['x'] += utils.add_noise(mean=constants.get_first_moment_err(),
                                  stdev=constants.get_first_moment_err_std(),
                                  shape=src['x'].shape,
                                  measurement=src['x'])
            src['y'] += utils.add_noise(mean=constants.get_first_moment_err(),
                                  stdev=constants.get_first_moment_err_std(),
                                  shape=src['y'].shape,
                                  measurement=src['y'])
            src['trace'] += utils.add_noise(mean=constants.get_second_moment_err(),
                                     stdev=constants.get_second_moment_err_std(),
                                     shape=src['trace'].shape,
                                     measurement=src['trace'])
        src['apMag'] = utils.flux_to_mag(src['modelFlux'], from_unit='nMgy')
        src['apMagErr'] = (2.5/np.log(10.0)) * src['apFluxErr'] / src['modelFlux']

        #####################################################
        # Final column renaming/reordering & saving to file #
        #####################################################
        src.rename(columns={'obsHistID': 'ccdVisitId',
                            'expMJD': 'MJD',
                            'FWHMeff': 'psf_fwhm',
                            'modelFlux': 'apFlux',
                            'mE1': 'e1',
                            'mE2': 'e2'}, inplace=True)
        src['e_final'], src['phi_final'] = utils.e1e2_to_ephi(src['e1'], src['e2'])
        src.drop(['mRrCc', 'offsetRa', 'offsetDec', 'fiveSigmaDepth'], axis=1, inplace=True)
        gc.collect()
        print("Number of observations: ", src['MJD'].nunique())
        print("Number of nonlenses: ", src['objectId'].nunique())

        src = src[self.source_columns]
        src.set_index('objectId', inplace=True)
        src.to_csv(save_file)
        end = time.time()

        print("Done making the source table with %d row(s) in %0.2f seconds using vectorization." %(len(src), end-start))

        self.sourceTable = src
        if self.DEBUG:
            return src


    #def make_source_table INHERITED
