import os, sys
# import galsim
import pandas as pd
import numpy as np

# PJM: if slrealizer package is installed, we should not need to
# adjust the sys.path here:
# realizer_path = os.path.join(os.environ['SLREALIZERDIR'], 'slrealizer')
# sys.path.insert(0, realizer_path)

# PJM: use slrealizer package:
# from realize_sdss import SDSSRealizer
# from utils.utils import *
from slrealizer import SDSSRealizer
from slrealizer.utils import *

if __name__=='__main__':
    """
    An annotated version of this script can be found in
    demo/Example+SLRealizer+Usage.ipynb.
    """
    # PJM: Avoid the need for an environment variable:
    # data_path = os.path.join(os.environ['SLREALIZERDIR'], 'data')
    data_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), '../data')

    catalog_f = os.path.join(data_path, 'sdss_processed.csv')
    observation_f = os.path.join(data_path, 'twinkles_observation_history.csv')

    output_nonlens_source_path = os.path.join(data_path, 'nonlens_source_table.csv')
    output_nonlens_object_path = os.path.join(data_path, 'nonlens_object_table.csv')

    db = pd.read_csv(catalog_f).sample(20, random_state=123).reset_index(drop=True)
    #db = pd.read_csv(catalog_f).sample(20000, random_state=123).reset_index(drop=True)
    # removing .query('mRrCc_u < 5.12')
    obs = pd.read_csv(observation_f)\
            .query("(expMJD < 65000) & (filter != 'y')")\
            .reset_index(drop=True)
    realizer = SDSSRealizer(observation=obs, catalog=db, debug=False, add_moment_noise=True, add_flux_noise=True)

    realizer.make_source_table_vectorized(save_path=output_nonlens_source_path)
    realizer.make_object_table(include_std=True,
                               source_table_path=output_nonlens_source_path,
                               object_table_path=output_nonlens_object_path)
