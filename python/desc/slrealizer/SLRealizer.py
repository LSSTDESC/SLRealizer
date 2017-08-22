#=====================================================
import numpy as np
import desc.slrealizer
import pandas as pd
import matplotlib.pyplot as plt
import pylab
import matplotlib
import math
import skimage
import random
import pandas
import corner
import dropbox
import extract_corner
import galsim_extract_corner
import om10
#from corner import corner
#=====================================================

class SLRealizer(object):

    """
    Contains the constructor and the key methods for SLRealizer module.
    Constructor reads in an OM10 catalog and observation history.
    Generates the toy catalog, plots the lensed system, and deblends sources using OM10 catalog and observation history.
    """

    def __init__(self, catalog=None, observation="../../../data/twinkles_observation_history.csv"):
        """
        Reads in a lens sample catalog and observation data.
        We assume lenses are OM10 lenses and observation file is .csv file
        """
        self.catalog = catalog
        self.observation = pd.read_csv(observation,index_col=0).as_matrix()
        filter = self.observation[:,1]
        self.observation_count = {'u': np.count_nonzero(filter=='u'),
                                      'g': np.count_nonzero(filter=='g'),
                                      'r': np.count_nonzero(filter=='r'),
                                      'i': np.count_nonzero(filter=='i'),
                                      'z': np.count_nonzero(filter=='z')}

    def plot_lens_random_date(self, lensID=None, convolve=False):
        """                                                                                                                      
                Given a specific lens, this code plots a lens after choosing a random observation epoch.                                           """
        if lensID is None:
            print 'No lens system selected for plotting.'
            return
        # Keep randomly selecting epochs until we get one that is not in the 'y' filter:
        filter = 'y'
        while filter == 'y':
            randomIndex = random.randint(0, 200)
            filter = self.observation[randomIndex][1]
        # Now visualize the lens system at the epoch defined by the randomIndex:
        desc.slrealizer.draw_model(self.observation[randomIndex],
                                   self.catalog.get_lens(lensID),
                                   convolve)
        return

    # after merging, change this one to deblend_test
    def deblend(self, lensID=None, null_deblend=True):
        """
        Given a lens system, this method deblends the source and plots the process of deblending.

        Parameters
        ---------
        lensID : int
        OM10 lens ID which can be used to identify the lensed system
        
        null_deblend : bool
        If true, assumes null deblender. Working deblender is currently not being supported
        """

        if lensID is None:
            print('No lens system selected for calculating the statistics')
            return
        if null_deblend is False:
            print('Sorry, working deblender is not being supported.')
            return
        # Keep randomly selecting epochs until we get one that is not in the 'y' filter:
        filter = 'y'
        while filter == 'y':
            randomIndex = random.randint(0, 200)
            filter = self.observation[randomIndex][1]
        image2 = desc.slrealizer.plot_all_objects(self.observation[randomIndex], self.catalog.get_lens(lensID))
        print('##################### PLOTTING ALL SOURCES ##################################')
        desc.slrealizer.show_color_map(image2)
        flux, first_moment_x, first_moment_y, covariance_matrix = desc.slrealizer.null_deblend(image2)
        print('first:', flux, first_moment_x, first_moment_y, covariance_matrix)
        image = desc.slrealizer.null_deblend_plot(flux, first_moment_x, first_moment_y, covariance_matrix)
        print('##################### AFTER NULL DEBLENDING ##################################')
        desc.slrealizer.show_color_map(image)

    def generate_cornerplot(self, color='black', object_table_dir = '../../../data/source_table.csv', option = None, params = None, range=None, galsim=False, overlap=None, normed=False, data=None, label=None):
        """
        Given a source table, this method plots the cornerplot. The user has to specify which attributes they want to look at.

        Parameters
        ----------
        option : list of string/strings
        Choose from size, position, color, ellipticity, and magnitude. The cornerplot will show the attribute that the user selected.
        If None, the user should specify the parameters(column names) in the source table.

        params : tuples of strings
        Names of columns in source_table that the user wants to see in the cornerplot.
        Only works when option is not specified.

        range : list of tuples for each plot
        Returns cornerplot (matplotlib.pyplot)
        """

        options = [None, 'size', 'x_position', 'y_position', 'color', 'ellipticity', 'magnitude', 'position', 'custom']
        object_table = pd.read_csv(object_table_dir)
        if ((option is None) and (params is None)):
            print('either specify params or option. You can choose among :')
            print(options)
            print('or specify columns in the source table that you want to see in the cornerplot.')
            return
        elif 'custom' in option:
            print('custom input')
        elif option is not None:
            data, label = np.array([]), []
            for elem in option:
                if elem not in options:
                    print(elem, 'is not in the option')
                    return
                method_name = 'calculate_'+elem
                cur_data, cur_label = getattr(galsim_extract_corner, method_name)(object_table)
                data = np.append(data, cur_data)
                label.extend(cur_label)
            data = data.reshape(len(label), len(object_table)).transpose()
        else:
            data, label = desc.slrealizer.extract_features(object_table, params)
        if overlap is None:
            fig = corner.corner(data, labels=label, color=color, smooth=1.0, range=range, hist_kwargs=dict(normed=normed))
        else:
            fig = corner.corner(data, labels=label, color=color, smooth=1.0, range=range, fig=overlap, hist_kwargs=dict(normed=normed))
        return fig

    def overlap(self, dir='../../../data/sdss.fits'):
        desc.slrealizer.save_as_catalog(catalog=dir)

    def make_source_catalog(self, dir='../../../data/source_catalog.csv', galsim=False):
        """
        Generates a full catalog(for each filter) of 200 lensed system and saves it 
        """
        self.catalog.select_random(maglim=23.3,area=20000.0,IQ=0.75)
        print('From the OM10 catalog, I am selecting LSST lenses')
        df = pd.DataFrame(columns=['MJD', 'filter', 'RA', 'RA_err', 'DEC', 'DEC_err', 'x', 'x_com_err', 'y', 'y_com_err', 'flux', 'flux_err', 'size', 'size_err', 'e', 'psf_sigma', 'sky', 'lensid'])
        ellipticity_upper_limit = desc.slrealizer.get_ellipticity_cut()
        debug_count = 0
        num_system = len(self.catalog.sample)
        print('number of system:', num_system)
        for j in xrange(200): # we will select 400 observation
            if self.observation[j][1] != 'y':
                for i in xrange(num_system): # we will use first 400 lenses
                    print('debug_count: ***************************** : ', debug_count)
                    debug_count += 1
                    if self.catalog.sample[i]['ELLIP'] < ellipticity_upper_limit: # ellipticity cut : 0.5
                        if galsim:
                            data = desc.slrealizer.galsim_generate_data(self.catalog.get_lens(self.catalog.sample[i]['LENSID']), self.observation[j])
                            if data is not None:
                                df.loc[len(df)]= data
                                print('data:', data)
                        else:
                            data = desc.slrealizer.generate_data(self.catalog.get_lens(self.catalog.sample[i]['LENSID']), self.observation[j])
                            df.loc[len(df)]= data
                            debug_count += 1
        df.set_index('lensid', inplace=True)
        df.to_csv(dir, index=True)
        if galsim:
            desc.slrealizer.dropbox_upload(dir, 'source_catalog_galsim.csv')
        else:
            desc.slrealizer.dropbox_upload(dir, 'source_catalog.csv')

    def make_object_catalog(self, source_table_dir='../../../data/source_catalog.csv', save_dir='../../../data/object_catalog.csv'):
        """
        From the source_table, make an object table by averaging the quantities for each filter and saves into the save_dir.
        """
        print('Reading in the catalog')
        df = pandas.read_csv(source_table_dir)
        lensID = df['lensid']
        lensID = lensID.drop_duplicates().as_matrix()
        column_name = ['lensid', 'g_flux', 'g_x', 'g_y', 'g_size', 'g_flux_err', 'g_x_com_err', 'g_y_com_err', 'g_size_err', 'g_e', 'z_flux', 'z_x', 'z_y', 'z_size', 'z_flux_err', 'z_x_com_err', 'z_y_com_err', 'z_size_err', 'z_e', 'i_flux', 'i_x', 'i_y', 'i_size', 'i_flux_err', 'i_x_com_err', 'i_y_com_err', 'i_size_err', 'i_e', 'r_flux', 'r_x', 'r_y', 'r_size', 'r_flux_err', 'r_x_com_err', 'r_y_com_err', 'r_size_err','r_e', 'u_flux', 'u_x', 'u_y', 'u_size', 'u_flux_err', 'u_x_com_err', 'u_y_com_err', 'u_size_err', 'u_e']
        source_table = pd.DataFrame(columns=column_name)
        for lens in lensID:
            lens_row = [lens]
            lens_array = df.loc[df['lensid'] == lens]
            print('**********LENSARRAY**********')
            print(lens_array)
            for filter in ['g', 'z', 'i', 'r', 'u']:
                print('*******INPUT****************')
                print(lens_array.loc[lens_array['filter'] == filter])
                print('*******RETURNED**************')
                print(desc.slrealizer.return_mean_properties(lens_array.loc[lens_array['filter'] == filter]))
                lens_row.extend(desc.slrealizer.return_mean_properties(lens_array.loc[lens_array['filter'] == filter]))
            source_table.loc[len(source_table)]= np.array(lens_row)
        source_table.to_csv(save_dir, index=False)
        desc.slrealizer.dropbox_upload(save_dir, 'object_catalog.csv')

    def test(self):
        desc.slrealizer.galsim_generate_data(self.catalog.get_lens(self.catalog.sample[0]['LENSID']), self.observation[0])
        print('*************!******************!*******************')
