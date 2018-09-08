"""
The :mod:`dataloader` module contains the :class:`Dataloader` utility class, which provides easy access to some standard test data files, including the OM10 lens database, and an example LSST observation history file. 
"""
from __future__ import absolute_import, division, print_function

import os

class Dataloader:
    """
    Utility class for reading in data files and doing some operations,
    such as querying or basic preprocessing, on them.

    Methods
    -------
    read(filename, is_test)
        Reads in the data file
    """

    def __init__(self):
        # Path of the data package where all data files reside
        from slrealizer import data
        self.data_dir = data.__path__[0]

    def read(self, filename, is_test):
        allowed_filenames = ['om10', 'observation', 'sdss']

        if filename=='om10':
            return self._read_om10(is_test=is_test)
        elif filename=='observation':
            return self._read_minion1060(is_test=is_test)
        elif filename=='sdss':
            return self._read_sdss(is_test=is_test)
        else:
            raise ValueError("Please choose one of %s" %(','.join(allowed_filenames)))

    def _read_om10(self, is_test):
        """Reads in the OM10 catalog.

        Reads in the OM10 catalog in the OM10's native database format, DB,
        and paints the lensed quasar systems in the catalog to have columns
        indicating multi-filter magnitudes

        Parameters
        ----------
        is_test : Boolean
            Whether to return a short test version of the catalog
            instead of the full catalog

        Returns
        -------
        DB
            The OM10 catalog.
        """
        from om10 import DB

        if is_test:
            catalog_path = os.path.join(self.data_dir, 'test_catalog.fits')
            lens_catalog = DB(catalog=catalog_path)
        else:
            lens_catalog = DB()

        lens_catalog.paint(synthetic=True)

        return lens_catalog

    def _read_minion1060(self, is_test):
        """Reads in the minion_1060 observational history.

        Reads in the minion_1060 observational history as a Pandas dataframe
        and queries out the Y filter.

        Parameters
        ----------
        is_test : Boolean

        Returns
        -------
        Pandas.Dataframe
            The minion_1060 catalog.
        """
        import pandas as pd

        minion1060_url = "https://www.dropbox.com/s/2fgjk6ip69d64kb/twinkles_observation_history.csv?dl=1"
        # Read in the catalog and query out the Y-filter, which does not exist in SDSS
        observation_catalog = pd.read_csv(minion1060_url).query("(filter != 'y')")

        if is_test:
            return observation_catalog.sample(20, random_state=123).reset_index(drop=True)
        else:
            return observation_catalog.reset_index(drop=True)

    def _read_sdss(self, is_test):
        """Reads in the SDSS object catalog.

        Reads in the processed SDSS object catalog as a Pandas dataframe.
        To see how the processing from the SDSS data release was done,
        see `slrealizer/Preprocessing+the+SDSS+Non-Lens+Catalog.ipynb`.

        Parameters
        ----------
        is_test : Boolean

        Returns
        -------
        Pandas.Dataframe
            The processed SDSS (non-lens) catalog.
        """
        import pandas as pd

        sdss_url = "https://www.dropbox.com/s/74moqzb5zhzmmiq/sdss_processed.csv?dl=1"
        sdss_catalog = pd.read_csv(sdss_url)

        if is_test:
            return sdss_catalog.sample(2, random_state=123).reset_index(drop=True)
        else:
            return sdss_catalog.reset_index(drop=True)
