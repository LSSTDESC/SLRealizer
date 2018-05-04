# ======================================================================
import os, unittest
import numpy as np
import scipy
import skimage
import desc.slrealizer
import pandas
import galsim
import matplotlib.pyplot as plt
import math
# ======================================================================

class TestCase(unittest.TestCase):
    
    """
    Tests the SLRealizer class.
    
    Notes
    -----
    Execute these tests with:
        nosetests
    from anywhere in the module, provided you have run
        pip install nose
    """
    
    # ------------------------------------------------------------------

    def setUp(self):
        # setup 2d array with given dimensions (declared as global variables)
        global x_min, x_max, y_min, y_max, distance
        x_min, x_max, y_min, y_max = desc.slrealizer.get_x_min(), desc.slrealizer.get_x_max(), desc.slrealizer.get_y_min(), desc.slrealizer.get_y_max()
        distance = desc.slrealizer.get_distance()
        x, y = np.mgrid[x_min:x_max:distance, y_min:y_max:distance]
        pos = np.dstack((x, y))
        number_of_rows = int((x_max-x_min)/distance)
        number_of_columns = int((y_max-y_min)/distance)
        image = [[0]*number_of_rows for _ in range(number_of_columns)]
        # return initialized 2d array
        return image

    def test_galsim(self):
        galaxy = galsim.Gaussian(flux=10.0, sigma=10)
        galaxy += galaxy
        """
        galaxy = galaxy.shift(1, 2)
        galaxy = galaxy.shear(e1=0.2)
        #psf = galsim.Gaussian(flux=1, sigma=3.0)
        #galaxy = galsim.Convolve([galaxy,psf]) 
        img = galaxy.drawImage(scale=.2)
        #shape_info = img.FindAdaptiveMom()
        first_moment_x, first_moment_y = (shape_info.moments_centroid.x-(len(img.array)/2.0))*0.2, (shape_info.moments_centroid.y-(len(img.array)/2.0))*0.2
        # total image intensity for best-fit elliptical Gaussian from adaptive moments. Normally, this field is simply equal to the image flux (for objects that follow a Gaussian light distribution, otherwise it is something approximating the flux). However, if the image was drawn using `drawImage(method='sb')` then moments_amp relates to the flux via flux=(moments_amp)*(pixel scale)^2.
        flux = shape_info.moments_amp
        # moment calculation needs to be changed, because 2d array is not returned                       
        I_xx = shape_info.moments_sigma # unit of pixels is returned, so change units for arcsec squared 
        print('I_xx', I_xx)
        print(shape_info.observed_shape.e1)
        print(shape_info.observed_shape.e2)
        print(shape_info.observed_shape.g1, shape_info.observed_shape.g2, shape_info.moments_n_iter)
        
        print('******************HEY****************')
        self.assertAlmostEqual(flux, 10)
        #self.assertAlmostEqual(first_moment_x, 1, places=0)
        #self.assertAlmostEqual(first_moment_y, 2, places=0)
        self.assertAlmostEqual(I_xx*0.2, 3.0, places=2)
        """
        img = galaxy.drawImage(scale=.2)   
        #random_number_generator = galsim.BaseDeviate(3) #random seed can be any number you want   
        #sky_level = math.pow(10, (22.5 - 23.5)/2.5)/5 # because Fb = 5 \sigma_b        
        #noise = galsim.PoissonNoise(random_number_generator, sky_level=sky_level) # I guess sky_level is just sky_sigma? 
        #print 'sky_level', sky_level, 'noise', noise
        #img.addNoise(noise)
        shape_info = img.FindAdaptiveMom()
        print('first moment', shape_info.moments_amp)
# ======================================================================

if __name__ == '__main__':

    unittest.main()

# ======================================================================
