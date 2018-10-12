from __future__ import absolute_import, division, print_function

from slrealizer import SLRealizer
import slrealizer.utils.utils as utils
import slrealizer.utils.moments as moments
import slrealizer.utils.constants as constants
import numpy as np
import pandas as pd
import galsim
import gc # need this to optimize memory usage

class CosmoDC2Realizer(SLRealizer):

    """
    A class that realizes objects in the CosmoDC2 catalog
    under the given observation conditions,
    into LSST-like "Object" and "Source" tables.
    """

    def __init__(self, observation, catalog, debug=False, add_moment_noise=True, add_flux_noise=True):
        #super(OM10Realizer, self).__init__(observation) # Didn't work for some reason
        self.as_super = super(CosmoDC2Realizer, self)
        self.as_super.__init__(observation, add_moment_noise=add_moment_noise, add_flux_noise=add_flux_noise)
        self.catalog = catalog
        self.num_systems = len(self.catalog.sample)
        self.DEBUG = debug
        
        self.input_centered_objects = None
        self.neighbors = None
        self.source_table = None
    
    def _preformat_input_objects(self, is_centered=1, catalog_to_preformat=None):
        """
        Preprocesses the input object catalog with the column conventions
        that can be interpreted by the class functions
        Parameters
        ==========
        is_centered: Bool
            whether the objects should be labeled as centered.
        new_object_catalog: Pandas.DataFrame
            if not given, defaults to the object catalog the class was instantiated with
        Returns
        =======
        preformatted_catalog: Pandas.DataFrame
            The resulting preformatted catalog.
        """
        # Set galaxy_id as index
        # Add disk_to_total_ratio_i = 1.0 - bulge_to_total_ratio_i
        # Add flux column for each filter
        # Add disk flux and bulge flux
        # Update magnitudes so that there are now disk magnitudes and bulge magnitudes
        # Separate df into bulge-related and disk-related
        # Add is_bulge column that takes value 1 if bulge and 0 if disk
        # Add is_centered = 1.0 (these objects have zero first moments by definition)
    
    def _preformat_source_table(self):
        """
        Initializes self.source_table by combining the observation and object catalogs
        with the column conventions that can be used by utility functions 
        related to moment estimation
        """
        # call _preformat_input_objects and assign output to self.input_centered_objects
        # optionally call _add_neighbors
        
        # Combine with observation catalog
        # Rename columns
        src.rename(columns={'obsHistID': 'ccdVisitId',
            'LENSID': 'objectId',
            'expMJD': 'MJD',
            'FWHMeff': 'psf_fwhm',
            'ELLIP': 'e',
            'PHIE': 'beta',}, inplace=True)
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
        
        self.source_table = src
    
    def _add_neighbors(self, neighbor_catalog, within_distance):
        """
        Adds neighbors from a separate catalog that are located 
        within a given radius, to self.input_centered_objects
        Parameters
        ==========
        neighbor_catalog: Pandas.DataFrame
            catalog containing neighbors of the centered objects
        within_radius: float
            radius in arcsec within which to include the neighboring objects
        """
        if self.input_centered_objects is None:
            raise ValueError("Must first define what objects will be at the center.")
        # Call _preformat_input_objects(is_centered=0, catalog_to_preformat=neighbor_catalog)
        # Combine group of centered + neighbors in a common key 

        
        