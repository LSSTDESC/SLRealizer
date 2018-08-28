from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import unittest
import os, sys
import shutil
import pandas as pd
import numpy as np

from slrealizer.realize_sdss import SDSSRealizer
from slrealizer.utils.utils import *

class SDSSRealizerTest(unittest.TestCase):

    """
    Tests the SDSSRealizer subclass.
    
    NOTE
    Execute these tests with:
        nosetests
    from anywhere in the module, provided you have run
        pip install nose
    """

    @classmethod
    def setUpClass(cls):
        from slrealizer import data
        
        # Input catalogs
        data_dir = data.__path__[0]
        #input_object_catalog = os.path.join(data_dir, 'sdss_processed.csv')
        #input_observation_catalog = os.path.join(data_dir, 'twinkles_observation_history.csv')
        input_object_catalog = "https://www.dropbox.com/s/74moqzb5zhzmmiq/sdss_processed.csv?dl=1"
        input_observation_catalog = "https://www.dropbox.com/s/2fgjk6ip69d64kb/twinkles_observation_history.csv?dl=1"

        # Output catalogs
        tests_dir = os.path.dirname(os.path.realpath(__file__))
        output_dir = os.path.join(tests_dir, 'test_output', 'test_sdssrealizer')
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.makedirs(output_dir)
        output_paths = {
        'rowbyrow_path': os.path.join(output_dir, 'rowbyrow_source.csv'),
        'vectorized_path': os.path.join(output_dir, 'vectorized_source.csv'),
        'object_path': os.path.join(output_dir, 'object.csv'),
        }

        for k, v in output_paths.items():
            setattr(cls, k, v)

        test_sdss_df = pd.read_csv(input_object_catalog).sample(2, random_state=123).reset_index(drop=True)
        test_obs_df = pd.read_csv(input_observation_catalog).query("(filter != 'y')").sample(20, random_state=123).reset_index(drop=True)
        # Instantiate SDSSRealizer
        cls.realizer = SDSSRealizer(observation=test_obs_df, catalog=test_sdss_df, debug=True, add_moment_noise=False, add_flux_noise=False)
        cls.nonlens_info = test_sdss_df.loc[0]
        cls.obs_info = test_obs_df.loc[0]

    def test_create_source_row(self):
        """ Tests whether create_source_row runs """
        self.realizer.create_source_row(lens_info=self.nonlens_info, obs_info=self.obs_info)
        
    def test_make_source_table(self):
        """ Tests whether make_source_table_rowbyrow runs """
        self.realizer.make_source_table_rowbyrow(save_file=self.rowbyrow_path)

    def test_make_source_table_vectorized(self):
        """ Tests whether make_source_table_vectorized runs """
        self.realizer.make_source_table_vectorized(save_file=self.vectorized_path)

    def test_make_object_table(self):
        """ Tests whether make_object_table runs """
        self.realizer.make_source_table_vectorized(save_file=self.vectorized_path)
        self.realizer.make_object_table(source_table_path=self.vectorized_path,
                                        object_table_path=self.object_path)

if __name__ == '__main__':
    unittest.main()
