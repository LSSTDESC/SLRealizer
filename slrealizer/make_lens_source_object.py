from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os, sys
from om10 import DB
import pandas as pd
import numpy as np
realizer_path = os.path.join(os.environ['SLREALIZERDIR'], 'slrealizer')
sys.path.insert(0, realizer_path)
from realize_om10 import OM10Realizer

if __name__=='__main__':
    """
    An annotated version of this script can be found in
    demo/Example+SLRealizer+Usage.ipynb.
    """
    data_path = os.path.join(os.environ['SLREALIZERDIR'], 'data')
    catalog_f = os.path.join(data_path, 'qso_mock.fits')
    observation_f = os.path.join(data_path, 'twinkles_observation_history.csv')

    output_lens_source_path = os.path.join(data_path, 'lens_source_table.csv')
    output_lens_object_path = os.path.join(data_path, 'lens_object_table.csv')

    db = DB(catalog=catalog_f)
    db.select_random(maglim=23.3, area=100.0, IQ=0.75)
    #db.select_random(maglim=23.3, area=1.e8, IQ=0.75)
    db.paint(synthetic=True)
    
    obs = pd.read_csv(observation_f)\
            .query("(expMJD < 65000) & (filter != 'y')")\
            .reset_index(drop=True)
    realizer = OM10Realizer(observation=obs, catalog=db, debug=False, add_moment_noise=True, add_flux_noise=True)

    realizer.make_source_table_vectorized(output_source_path=output_lens_source_path,
                                          include_time_variability=True)
    
    realizer.make_object_table(include_std=True,
                               source_table_path=output_lens_source_path,
                               object_table_path=output_lens_object_path)
