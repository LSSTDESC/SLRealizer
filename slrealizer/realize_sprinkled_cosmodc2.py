from __future__ import absolute_import, division, print_function

from slrealizer import SLRealizer
import slrealizer.utils.utils as utils
import slrealizer.utils.moments as moments
import slrealizer.utils.constants as constants
import numpy as np
import pandas as pd
import galsim
import gc # need this to optimize memory usage

class SprinkledCosmoDC2Realizer(CosmoDC2Realizer):

    """
    A class that sprinkles mock lenses in the CosmoDC2 catalog
    and realizes them under the given observation conditions,
    into LSST-like "Object" and "Source" tables.
    """

    def __init__(self, observation, catalog, debug=False, add_moment_noise=True, add_flux_noise=True):
        #super(OM10Realizer, self).__init__(observation) # Didn't work for some reason
        self.as_super = super(SprinkledCosmoDC2Realizer, self)
        self.as_super.__init__(observation, add_moment_noise=add_moment_noise, add_flux_noise=add_flux_noise)
        self.catalog = catalog
        self.num_systems = len(self.catalog.sample)
        self.DEBUG = debug
        
        self.input_centered_objects = None
        self.neighbors = None
        self.source_table = None
    
    # INHERITS def _preformat_input_objects
    # INHERITS def _preformat_source_table
    # INHERITS _add_neighbors
    
    def assign_stellar_mass():
        """
        Assigns stellar mass to OM10 objects based on their velocity dispersion
        """
        pass
        
    def assign_sed_radius():
        """
        Matches OM10 objects with those in the extragalactic catalog with
        similar stellar mass, redshift, and ellipticity and gets the
        SED and radius
        """
        pass
        
    def cache_agn_to_replace():
        """
        Goes through all the AGNs in the extragalactic catalog and caches the IDs 
        that will be replaced
        """
        pass
    
    def sprinkle_lenses():