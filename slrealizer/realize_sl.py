from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
from slrealizer.utils.utils import *
from slrealizer.utils.constants import *
import pandas as pd
import random
import galsim

class SLRealizer(object):

    """
   
    Class equipped with utility functions
    for making the LSST-like object and source tables
    inherited by child classes which are associated with a specific non-LSST catalog, 
    e.g. the child class OM10Realizer converts the OM10 catalog into a mock LSST catalog
    
    """

    def __init__(self, observation, add_moment_noise, add_flux_noise):
        """
        Reads in a lens sample catalog and observation data.
        We assume lenses are OM10 lenses and observation file is a pandas df
        """
        self.observation = observation
        self.num_obs = len(self.observation)
        
        # GalSim drawImage params
        self.fft_params = galsim.GSParams(maximum_fft_size=10240)
        self.pixel_scale = 0.1
        self.nx, self.ny = 49, 49 
        
        # Source table df
        self.source_table = None
        # Source table column list
        self.source_columns = ['MJD', 'ccdVisitId', 'objectId', 'filter', 'psf_fwhm', 'x', 'y', 'apFlux', 'apFluxErr', 'apMag', 'apMagErr', 'trace', 'e1', 'e2', 'e_final', 'phi_final', ]
        
        # Controlling randomness
        self.add_moment_noise = add_moment_noise
        self.add_flux_noise = add_flux_noise
        if not self.add_flux_noise and not self.add_moment_noise:
            self.remove_random = True
        else:
            self.remove_random = False
        self.seed = 123
        np.random.seed(self.seed)
        
    def get_obs_info(self, obsID=None, rownum=None):
        if obsID is not None and rownum is not None:
            raise ValueError("Need to define either obsID or rownum, not both.")
        
        if obsID is not None:
            return self.observation.loc[self.observation['obsHistID']==obsID]
        elif rownum is not None:
            return self.observation.loc[rownum]
    
    def _get_lens_info(self, lensID):
        ''' This function will depend on the format of each lens catalog '''
        raise NotImplementedError
        
    def draw_system(self, galsimInput, obs_info, save_path=None):
        '''
        Draws all objects of the given lens system
        in the given observation conditions using GalSim

        Keyword arguments:
        lens_info -- a row of the OM10 DB
        obs_info -- a row of the observation history df
        save_path -- path in which to save the image

        Returns:
        A GalSim object of the aggregate system used to render
        the image
        '''
        histID, MJD, band, PSF_FWHM, sky_mag = obs_info

        # Lens galaxy
        galaxy = galsim.Gaussian(half_light_radius=galsimInput['half_light_radius'],\
                                 flux=galsimInput['flux'])\
                       .shear(e=galsimInput['e'], beta=galsimInput['beta'])
        # Lensed quasar
        for i in xrange(galsimInput['num_objects']):
            lens = galsim.Gaussian(flux=galsimInput['flux_'+str(i)], sigma=0.0)\
                         .shift(galsimInput['xy_'+str(i)])
            galaxy += lens
            
        psf = galsim.Gaussian(flux=1.0, fwhm=PSF_FWHM)
        galsim_obj = galsim.Convolve([galaxy, psf], gsparams=self.fft_params)
        galsim_img = galsim_obj.drawImage(nx=self.nx, ny=self.ny, scale=self.pixel_scale)
        if save_path is not None:
            plt.imshow(galsim_img.array, interpolation='none', aspect='auto')
            plt.savefig(save_path)
            plt.close()
        return galsim_img
        
    def estimate_parameters(self, galsim_img, method="raw_numerical"):
        """
        Performs shape estimati on on the galsim_img 
        using either GalSim's HSM shape estimator or 
        a native numerical moment calculator
        under the observation conditions obs_info
        
        Keyword arguments:
        galsim_img -- GalSim's Image object on which parameters will be estimated
        method -- one of "hsm" (GalSim's HSM shape estimator) or 
                  "raw_numerical" (a native numerical moment calculator) [default: "raw"]
        
        Returns
        a dictionary of the lens properties, 
        which can be used to draw the emulated image
        """
        estimated_params = {}
        if method == "hsm":       
            try:
                shape_info = galsim_img.FindAdaptiveMom(guess_sig=self.pixel_scale*10.0)
            except:
                #print "HSM failed"
                return None

            # Calculate the real position from the arbitrary pixel position
            pixelCenter = galsim.PositionD(x=shape_info.moments_centroid.x, y=shape_info.moments_centroid.y)
            estimated_params['x'], estimated_params['y'] = pixel_to_physical(shape_info.moments_centroid.x, self.nx, self.pixel_scale),\
                                             pixel_to_physical(shape_info.moments_centroid.y, self.ny, self.pixel_scale)

            estimated_params['apFlux'] = float(np.sum(galsim_img.array))
            if self.DEBUG:
                estimated_params['hlr'] = galsim_img.calculateHLR(center=pixelCenter)
                estimated_params['det'] = (shape_info.moments_sigma * self.pixel_scale)**4.0
            estimated_params['trace'] = 2.0*galsim_img.calculateMomentRadius(center=pixelCenter, rtype='trace')**2.0
            estimated_params['e1'] = shape_info.observed_shape.e1
            estimated_params['e2'] = shape_info.observed_shape.e2
            estimated_params['e_final'] = shape_info.observed_shape.e
            estimated_params['phi_final'] = shape_info.observed_shape.beta
        elif method == "raw_numerical":
            image_array = galsim_img.array
            Ix, Iy = get_first_moments_from_image(image_array, self.pixel_scale)
            Ixx, Ixy, Iyy = get_second_moments_from_image(image_array, self.pixel_scale)
            estimated_params['apFlux'] = np.sum(image_array)
            estimated_params['x'] = Ix
            estimated_params['y'] = Iy
            trace = Ixx + Iyy
            estimated_params['e1'] = (Ixx - Iyy)/trace
            estimated_params['e2'] = 2.0*Ixy/trace
            estimated_params['trace'] = trace
            estimated_params['e_final'], estimated_params['phi_final'] = e1e2_to_ephi(estimated_params['e1'], estimated_params['e2'])
        else:
            raise ValueError("Please enter a valid method, either 'hsm' or 'raw_numerical'")

        return estimated_params
    
    def draw_emulated_system(self, estimated_params):
        """
        Draws the emulated system, i.e. draws the aggregate system
        from properties HSM derived from the image, which was in turn
        drawn from the catalog's truth properties.
        Only runs when DEBUG == True.
        
        Returns
        a GalSim Image object of the emulated system
        """
        if not self.DEBUG:
            raise ValueError("Only runs in debug mode")
        system = galsim.Gaussian(flux=estimated_params['apFlux'], half_light_radius=estimated_params['hlr'])\
                       .shift(float(estimated_params['x']), float(estimated_params['y']))\
                       .shear(e1=estimated_params['e1'], e2=estimated_params['e2'])
        emulatedImg = system.drawImage(nx=self.nx, ny=self.ny, scale=self.pixel_scale, method='no_pixel')
        return emulatedImg
    
    def create_source_row(self, derived_params, objectId, obs_info):
        '''
        Returns a dictionary of lens system's properties
        computed the image of one lens system and the observation conditions,
        which makes up a row of the source table.

        Keyword arguments:
        derived_params -- derived lens properties 
        obs_info -- a row of the observation history df
        
        Returns
        A dictionary with properties derived from HSM estimation
        (See code for which properties)
        '''
        histID, MJD, band, PSF_FWHM, sky_mag = obs_info
        
        derived_params['apFluxErr'] = mag_to_flux(sky_mag-22.5)/5.0 # because Fb = 5 \sigma_b
        if self.add_moment_noise:
            derived_params['trace'] += add_noise(mean=get_second_moment_err(), 
                                               stdev=get_second_moment_err_std(), 
                                               measurement=derived_params['trace'])
            derived_params['x'] += add_noise(mean=get_first_moment_err(), 
                                           stdev=get_first_moment_err_std(), 
                                           measurement=derived_params['x'])
            derived_params['y'] += add_noise(mean=get_first_moment_err(), 
                                           stdev=get_first_moment_err_std(), 
                                           measurement=derived_params['y']) 
        if self.add_flux_noise:
            derived_params['apFlux'] += add_noise(mean=0.0, 
                                                stdev=derived_params['apFluxErr']) # flux rms not skyErr
        derived_params['apMag'] = flux_to_mag(derived_params['apFlux'], from_unit='nMgy')
        derived_params['apMagErr'] = (2.5/np.log(10.0)) * derived_params['apFluxErr']/derived_params['apFlux']
        
        row = {'MJD': MJD, 'ccdVisitId': histID, 'filter': band, 'x': derived_params['x'], 'y': derived_params['y'],
               'apFlux': derived_params['apFlux'], 'apFluxErr': derived_params['apFluxErr'], 
               'apMag': derived_params['apMag'], 'apMagErr': derived_params['apMagErr'],
               'trace': derived_params['trace'],
               'e1': derived_params['e1'], 'e2': derived_params['e2'], 
               'e_final': derived_params['e_final'], 'phi_final': derived_params['phi_final'],
               'psf_fwhm': PSF_FWHM, 'objectId': objectId}
        return row

    def make_source_table_rowbyrow(self, save_file, method="analytical"):
        import time
        """
        Returns a source table generated from all the lens systems in the catalog
        under all the observation conditions in the observation history,
        and saves it as a .csv file.

        Keyword arguments:
        save_file -- path into which output source table will be saved
        method -- how to calculate moments for each row
                  (See method estimate_parameters for details about each option)         
        """
        start = time.time()
        print("Began making the source catalog.")
        
        df = pd.DataFrame(columns=self.source_columns)
        #ellipticity_upper_limit = desc.slrealizer.get_ellipticity_cut()
        print("Number of systems: %d, number of observations: %d" %(self.num_systems, self.num_obs))
        
        hsm_failed = 0
        
        for j in xrange(self.num_obs):
            for i in xrange(self.num_systems):
                row = self.create_source_row(lens_info=self.get_lens_info(rownum=i),
                                             obs_info=self.observation.loc[j],
                                             method=method)
                if row == None:
                    hsm_failed += 1
                else:
                    df = df.append(row, ignore_index=True)
        
        df = df.infer_objects()
        df = df[self.source_columns]
        df.set_index('objectId', inplace=True)
        df.to_csv(save_file, index=True)
        
        end = time.time()
        if method == 'hsm':
            print("Done making the source table which has %d row(s) in %0.2f hours, after getting %d errors from HSM failure." %(len(df), (end - start)/3600.0, hsm_failed))
        else:
            print("Done making the source table with %s method in %0.2f minutes." %(method, (end - start)/60.0))
#        desc.slrealizer.dropbox_upload(dir, 'source_catalog_new.csv')

        self.sourceTable = df
        if self.DEBUG:
            return df

    def make_object_table(self, object_table_path, source_table_path=None, include_std=False):

        """
        Generates the object table from the given source table at source_table_path
        by averaging the properties for each filter, and saves it as object_table_path.
        """
        import time
        import gc
        
        if object_table_path is None:
            raise ValueError("Must provide save path of the output object table.")
        
        if source_table_path is not None:
            print("Reading in the source table at %s ..." %source_table_path)
            obj = pd.read_csv(source_table_path)
            obj.set_index('objectId', inplace=True)
        elif self.sourceTable is not None:            
            print("Reading in Pandas Dataframe of most recent source table generated... ")
            obj = self.sourceTable.copy()
        else:
            raise ValueError("Must provide a source table path or generate a source table at least once using this Realizer object.")

        start = time.time()
        
        obj.drop(['ccdVisitId', 'psf_fwhm'], axis=1, inplace=True)
        # Define (filter-nonspecific) properties to go in object table columns
        keepCols = list(obj.columns.values)
        keepCols.remove('filter')
        # Pivot the filter values to the column axis
        obj = obj.pivot_table(values=keepCols, index=[obj.index, 'MJD'], columns='filter', aggfunc='first')
        # Collapse multi-indexed column using filter_property formatting
        obj.columns = obj.columns.map('{0[1]}_{0[0]}'.format)
        gc.collect()
        
        #if self.DEBUG: print(obj.columns.values)
        
        # Take mean, optional std of properties across observed times for each object
        obj.reset_index(inplace=True)
        obj.drop('MJD', axis=1, inplace=True)
        obj = obj.groupby('objectId', sort=False)
        means = obj.mean()
        if include_std:
            stds = obj.std()
            obj = means.join(stds, lsuffix='', rsuffix='-std')
        else:
            obj = means
        gc.collect()
        
        # Drop examples with missing values
        obj.dropna(how='any', inplace=True)
        # Get x, y values relative to the r-band
        for b in 'ugriz':
            obj[b + '_' + 'x'] = obj[b + '_' + 'x'] - obj['r_x']
            obj[b + '_' + 'y'] = obj[b + '_' + 'y'] - obj['r_y']
        end = time.time()
        
        # Save as csv file
        obj.to_csv(object_table_path, index=False)
        print("Done making the object table in %0.2f seconds." %(end-start))
        #if self.DEBUG:
            #print("Object table columns: ", obj.columns)

        #desc.slrealizer.dropbox_upload(save_dir, 'object_catalog_new.csv') #this uploads to the desc account
    
    def include_quasar_variability(self, save_output=False, input_source_path=None, output_source_path=None):
        """
        Takes a source table and adds the intrinsic variability of the quasar images
        using the generative model introduced in MacLeod et al (2010)
        
        Keyword arguments:
        save_output -- whether to save the output to disk [default: False]
        input_source_path -- path of input source table to be altered [default: None]
        output_source_path -- path of output source table containing time variability [default: None]
        """
        
        import gc
        import time
        
        start = time.time()
        if input_source_path is None:
            try:
                print("Reading in the most recent source table...")
                src = self.source_table
            except ValueError:
                print("Realizer has not generated a source table yet.")
        else:
            try:
                print("Reading in the source table at %s" %input_source_path)
                src = pd.read_csv(input_source_path)
            except ValueError:
                print("Please input a valid path to the source table.")
            
        if self.DEBUG:
            NUM_TIMES = src['MJD'].nunique()
            NUM_OBJECTS = src['objectId'].nunique()
            print(NUM_TIMES, NUM_OBJECTS)
        
        for q in range(4):
            magnitude_type = 'q_mag_' + str(q)
            # Filter-specific dictionary of apMag values to replace src apMag
            apMag_filter = {}
            for b in 'ugriz':
                #########################
                # Adding useful columns #
                #########################
                # Filter-specific version of src with only the columns we need
                src_timevar = src.query('filter == @b')[['filter', 'objectId', 'MJD', magnitude_type, ]].copy()
                # Total number of time steps in this filter
                NUM_TIMES_FILTER = src_timevar['MJD'].nunique()
                # Set MJD relative to zero
                src_timevar['MJD_diff'] = src_timevar['MJD'] - np.min(src_timevar['MJD'].values)
                # Map ccdVisitId to integers starting from 0
                sorted_obs_id = np.sort(src_timevar['MJD'].unique())
                time_map = dict(zip(sorted_obs_id, range(NUM_TIMES_FILTER)))
                src_timevar['time_index'] = src_timevar['MJD'].map(time_map).astype(int)
                # Add a column of time elapsed since last observation, d_time
                src_timevar.sort_values(['objectId', 'filter', 'MJD', 'MJD_diff'], axis=0, inplace=True)
                src_timevar['d_time'] = src_timevar['MJD_diff'] - src_timevar['MJD_diff'].shift(+1)
                src_timevar['d_time'].fillna(0.0, inplace=True)
                src_timevar['d_time'] = np.clip(src_timevar['d_time'], a_min=0.0, a_max=None)
                gc.collect()

                ##############################
                # Computing time variability #
                ##############################
                # Parameters of the generative model (hand-picked)
                MU = 0.0
                TAU = 20.0 #np.power(10.0, 2.4) # days
                S_INF = 0.14 # mag
                # Add a column to store the variability and initialize to zero
                src_timevar['intrinsic_mag'] = 0.0
                # Pivot to get time sequence to be horizontal
                src_timepivot = src_timevar.pivot_table(index=['filter', 'objectId'], 
                                                        columns=['time_index',], 
                                                        values=['d_time', magnitude_type, 'intrinsic_mag', ])
                # Update 'intrinsic_mag' column
                for t in range(1, NUM_TIMES_FILTER):
                    src_timepivot.loc[b, :]['intrinsic_mag'][t] = np.random.normal(loc=src_timepivot['intrinsic_mag'][t - 1].values*np.exp(-src_timepivot['d_time'][t].values/TAU) + MU*(1.0 - np.exp(-src_timepivot['d_time'][t].values/TAU)),
                                                                                   scale=0.5*S_INF**2.0*(1.0 - np.exp(-2.0*src_timepivot['d_time'][t].values/TAU)))

                ########################
                # Final column sorting #
                ########################
                # Add computed variability to magnitude_type column
                src_timepivot[magnitude_type] = src_timepivot[magnitude_type] + src_timepivot['intrinsic_mag']    
                # Pivot the time sequence back, to be vertical
                src_timepivot = src_timepivot.stack(dropna=False)
                # Drop columns we'll no longer use
                src_timepivot.drop(['intrinsic_mag', 'd_time'], axis=1, inplace=True)
                # As precaution, sort before storing updated apMag
                src_timepivot.reset_index(inplace=True)
                src_timepivot.sort_values(by=['filter', 'objectId', 'time_index'], inplace=True)
                # Store filter-specific apMag values 
                apMag_filter[b] = src_timepivot[magnitude_type].values
                gc.collect()

            # Update apMag in source table, filter by filter
            for b in 'ugriz':
                src.loc[src['filter']==b, magnitude_type] = apMag_filter[b]
            
        gc.collect()
        if self.DEBUG:
            print("Result of adding time variability: ")
            print("Number of observations: ", src['MJD'].nunique())
            print("Number of objects: ", src['objectId'].nunique())
        
        src.set_index('objectId', inplace=True)
        end = time.time()
        
        print("Done adding time variability with %d row(s) in %0.2f seconds using vectorization." %(len(src), end-start))
        if save_output:
            print("Saving the new source table with time variability at %s" %output_source_path)
            src.to_csv(output_source_path)
            
        self.source_table = src
    
    def _include_moments(self, inplace=True, input_dict=None):
        """
        Adds columns of first and second moments (analytically computed)
        to self.source_table (if inplace=True) or some other dictionary

        Keyword arguments:
        inplace -- whether to alter the self.source_table. 
                   If input_dict is None, will be assumed to be True. [default: True]
        input_dict -- object of type dict that will be used to compute the moments.
                      If inplace is True, will be ignored. [default: None]

        Returns (only if inplace is False):
        a new dictionary that is a version of the original input_dict with
        moment-related columns added
        """
        
        return_dict = False
        if inplace or input_dict is None:
            src = self.source_table
        elif input_dict is not None:
            is_dictionary = (type(input_dict) is dict)
            valid_columns = ['lens_flux', 'apFlux', 'e', 'beta', 'psf_fwhm'] 
            valid_columns += [value + '_%d' %idx for value in ['q_flux', 'XIMG', 'YIMG'] for idx in range(4)]
            has_valid_columns = all(col in input_dict for col in valid_columns)
            if is_dictionary and has_valid_columns:
               src = input_dict
               return_dict = True
            else:
                raise ValueError("Keyword input_dict must be a dictionary and contain all the required keys.")
        
        # Calculate flux ratios (for weighted moments)
        src['lensFluxRatio'] = src['lens_flux']/src['apFlux']
        for q in range(4):
            src['qFluxRatio_' + str(q)] = src['q_flux_' + str(q)]/src['apFlux']
        if not return_dict:
            src.drop(['lens_flux'] + ['q_flux_' + str(q) for q in range(4)], axis=1, inplace=True)
        
        #################
        # FIRST MOMENTS #
        #################
        src['x'], src['y'] = 0.0, 0.0
        for q in range(4):
            src['x'] += src['qFluxRatio_' + str(q)]*src['XIMG_' + str(q)]
            src['y'] += src['qFluxRatio_' + str(q)]*src['YIMG_' + str(q)]
        if self.add_moment_noise:
            src['x'] += add_noise(mean=get_first_moment_err(), 
                                  stdev=get_first_moment_err_std(), 
                                  shape=src['x'].shape,
                                  measurement=src['x'])
            src['y'] += add_noise(mean=get_first_moment_err(), 
                                  stdev=get_first_moment_err_std(), 
                                  shape=src['y'].shape,
                                  measurement=src['y'])
        ##################
        # SECOND MOMENTS #
        ##################
        # Read in lens shape parameters
        src['minor_to_major'] = np.power((1.0 - src['e'])/(1.0 + src['e']), 0.5) # q parameter in galsim.shear
        src['beta'] = np.radians(src['beta']) # beta parameter in galsim.shear 
         # Arbitrarily set REFF_T to 1.0
        #src['sigmasq_lens'] = np.power(hlr_to_sigma(src['REFF_T']), 2.0)
        src['sigmasq_lens'] = np.power(hlr_to_sigma(1.0), 2.0)
        
        # Initialize with lens contributions
        src['lam1'] = src['sigmasq_lens']/src['minor_to_major']
        src['lam2'] = src['sigmasq_lens']*src['minor_to_major']
        src['lens_Ixx'] = src['lam1']*np.power(np.cos(src['beta']), 2.0) + src['lam2']*np.power(np.sin(src['beta']), 2.0)
        src['lens_Iyy'] = src['lam1']*np.power(np.sin(src['beta']), 2.0) + src['lam2']*np.power(np.cos(src['beta']), 2.0)
        src['lens_Ixy'] = (src['lam1'] - src['lam2'])*np.cos(src['beta'])*np.sin(src['beta'])
        src['Ixx'] = src['lensFluxRatio']*(src['lens_Ixx'] + np.power(src['x'], 2.0)) 
        src['Iyy'] = src['lensFluxRatio']*(src['lens_Iyy'] + np.power(src['y'], 2.0))
        src['Ixy'] = src['lensFluxRatio']*(src['lens_Ixy'] - src['x']*src['y'])
        if not return_dict:
            src.drop(['lam1', 'lam2', 'lens_Ixx', 'lens_Iyy', 'lens_Ixy'], axis=1, inplace=True)
        # Add quasar contributions
        for q in range(4):
            src['Ixx'] += src['qFluxRatio_' + str(q)]*(np.power(src['XIMG_' + str(q)] - src['x'], 2.0))
            src['Iyy'] += src['qFluxRatio_' + str(q)]*(np.power(src['YIMG_' + str(q)] - src['y'], 2.0))
            src['Ixy'] += src['qFluxRatio_' + str(q)]*(src['XIMG_' + str(q)] - src['x'])\
                                                     *(src['YIMG_' + str(q)] - src['y'])
        # Add PSF
        src['sigmasq_psf'] = np.power(fwhm_to_sigma(src['psf_fwhm']), 2.0)
        src['Ixx'] += src['sigmasq_psf']
        src['Iyy'] += src['sigmasq_psf']
                
        # Get trace and ellipticities
        src['trace'] = src['Ixx'] + src['Iyy']
        if self.add_moment_noise:
            src['trace'] += add_noise(mean=get_second_moment_err(), 
                                      stdev=get_second_moment_err_std(), 
                                      shape=src['trace'].shape,
                                      measurement=src['trace'])
        src['e1'] = (src['Ixx'] - src['Iyy'])/src['trace']
        src['e2'] = 2.0*src['Ixy']/src['trace']
        src['e_final'], src['phi_final'] = e1e2_to_ephi(src['e1'], src['e2'])
        
        if inplace:
            self.source_table = src

        if return_dict:
            return src
    
    def compare_truth_vs_emulated(self, lensID=None, rownum=None, save_dir=None):
        """                                                                                                                   
        Draws two images of the lens system with the given rownum
        under randomly chosen observation conditions,
        one from the catalog info (truth image) and
        another from HSM's estimation of the truth image

        Keyword arguments:
        - lensID: an integer specifying the ID of the lens system in the OM10 catalog

        Returns:
        A tuple of
        - the truth image drawn from the catalog info
        - the emulated image drawn from HSM's derived info of aggregate system 

        """
        # Only works in debug mode.
        self.debug = True
        # Randomly select observation ID
        obs_rownum = random.randint(0, self.num_obs)
        
        # Render the truth image
        truth_img = self.draw_system(self.get_lens_info(lensID=lensID, rownum=rownum),\
                                    self.get_observation(rownum=obs_rownum),\
                                    save_dir)
        
        # Render the emulated image under observation conditions indexed by obs_rownum
        fig, axes = plt.subplots(2, figsize=(5, 10))
        axes[0].imshow(truth_img.array, interpolation='none', aspect='auto')
        axes[0].set_title("TRUE MODEL IMAGE")
        estimated_params = self.estimate_hsm(image=truth_img, observation=self.observation.loc[obs_rownum])
        emulated_img = self.draw_emulated_system(estimated_params)
        axes[1].imshow(img.array, interpolation='none', aspect='auto')
        axes[1].set_title("EMULATED IMAGE")
        if save_dir is not None:
            plt.savefig(save_dir + 'truth_vs_emulated.png')
        return truth_img, emulated_img
