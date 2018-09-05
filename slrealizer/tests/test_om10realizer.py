from __future__ import absolute_import, division, print_function

import unittest
import os, sys
import shutil
import pandas as pd
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from om10 import DB
from slrealizer import OM10Realizer

class OM10RealizerTest(unittest.TestCase):

    """
    Tests the OM10Realizer subclass.
    
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
        input_object_catalog = os.path.join(data_dir, 'test_catalog.fits')
        #input_observation_catalog = os.path.join(data_dir, 'twinkles_observation_history.csv')
        input_observation_catalog = "https://www.dropbox.com/s/2fgjk6ip69d64kb/twinkles_observation_history.csv?dl=1"

        # Output catalogs
        tests_dir = os.path.dirname(os.path.realpath(__file__))
        output_dir = os.path.join(tests_dir, 'test_output', 'test_om10realizer')
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.makedirs(output_dir)
        output_paths = {
        'rowbyrow_analytical_path': os.path.join(output_dir, 'rowbyrow_ana_source.csv'),
        'rowbyrow_hsm_numerical_path': os.path.join(output_dir, 'rowbyrow_hsm_num_source.csv'),
        'rowbyrow_raw_numerical_path': os.path.join(output_dir, 'rowbyrow_raw_num_source.csv'),
        'vectorized_path': os.path.join(output_dir, 'vectorized_source.csv'),
        'object_path': os.path.join(output_dir, 'object.csv'),
        }

        for k, v in output_paths.items():
            setattr(cls, k, v)

        test_om10_db = DB(catalog=input_object_catalog)
        test_om10_db.paint(synthetic=True)
        test_obs_df = pd.read_csv(input_observation_catalog).query("(filter != 'y')").sample(20, random_state=123).reset_index(drop=True)
        # Instantiate OM10Realizer
        cls.realizer = OM10Realizer(observation=test_obs_df, catalog=test_om10_db, debug=True, add_moment_noise=False, add_flux_noise=False) 
        cls.lens_info = test_om10_db.sample[0]
        cls.obs_info = test_obs_df.loc[0]

    def test_om10_to_galsim(self):
        """ Tests whether _om10_to_galsim method runs """
        self.realizer._om10_to_galsim(lens_info=self.lens_info, band=self.obs_info['filter'])

    def test_draw_system(self):
        """ Tests whether draw_system method runs """ 
        self.realizer.draw_system(lens_info=self.lens_info, obs_info=self.obs_info)

    def test_estimate_hsm(self):
        """ Tests whether estimate_parameters method runs """ 
        self.realizer.estimate_parameters(lens_info=self.lens_info, obs_info=self.obs_info)

    def test_om10_to_lsst(self):
        """ Tests whether _om10_to_lsst runs """
        self.realizer._om10_to_lsst(lens_info=self.lens_info, obs_info=self.obs_info)

    def test_create_source_row(self):
        """ 
        Tests whether create_source_row runs with each of the options 
        for moment calculation method 
        """
        for m in ["analytical", "hsm", "raw_numerical"]:
            self.realizer.create_source_row(lens_info=self.lens_info, obs_info=self.obs_info, method=m)

    def test_make_source_table_rowbyrow_numerical(self):
        """ 
        Tests whether make_source_table_rowbyrow runs with two numerical options 
        for moment calculation method, 'hsm' and 'raw_numerical'
        """
        self.realizer.make_source_table_rowbyrow(save_file=self.rowbyrow_raw_numerical_path, method="raw_numerical")
        self.realizer.make_source_table_rowbyrow(save_file=self.rowbyrow_hsm_numerical_path, method="hsm")

    def test_make_source_table_analytical(self):
        """ 
        Tests whether make_source_table_vectorized run and
        whether make_source_table_rowbyrow runs with method = 'analytical',
        and checks whether the two output source tables have the same values
        """
        rowbyrow = self.realizer.make_source_table_rowbyrow(save_file=self.rowbyrow_analytical_path, method="analytical")
        vectorized = self.realizer.make_source_table_vectorized(output_source_path=self.vectorized_path, include_time_variability=False)

        #print(rowbyrow.dtypes, vectorized.dtypes)

        # Select only the columns with numeric values, and turn into Numpy array
        rowbyrow_float = rowbyrow.select_dtypes(include=[np.number]).values
        vectorized_float = vectorized.select_dtypes(include=[np.number]).values

        self.assertTrue(np.allclose(rowbyrow_float, vectorized_float, rtol=1e-05, atol=1e-05))

    def test_make_object_table(self):
        """ Tests whether make_object_table runs """
        self.realizer.make_source_table_vectorized(output_source_path=self.vectorized_path, include_time_variability=False)
        self.realizer.make_object_table(include_std=False, source_table_path=self.vectorized_path, object_table_path=self.object_path)

if __name__ == '__main__':
    unittest.main()
