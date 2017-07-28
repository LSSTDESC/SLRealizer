# ======================================================================                                    
from __future__ import print_function
from scipy.stats import multivariate_normal
#import plotting
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, MaxNLocator
from numpy import linspace
import matplotlib
import scipy
import scipy.optimize as opt
import math
import photutils
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils import CircularAperture
from photutils import DAOStarFinder
import moment
import desc.slrealizer
import scipy.stats
from scipy.stats import moment
from scipy import ndimage
import skimage
import skimage.measure
from skimage.measure import moments

# ======================================================================                                    
#https://github.com/scikit-image/scikit-image/blob/master/skimage/measure/_moments_cy.pyx


#global variable that controls the size of the plot
x_min = -5.0
x_max = 5.0
y_min = -5.0
y_max = 5.0
distance = 0.01

number_of_rows = int((x_max - x_min)/distance)
number_of_columns = int((y_max - y_min)/distance)

def null_deblend_v3(image2):
    x, y = np.mgrid[x_min:x_max:distance, y_min:y_max:distance]
    pos = np.dstack((x, y))
    #moment_matrix = skimage.measure.moments(image2)
    zeroth_moment = desc.slrealizer.zeroth_moment(image2)
    first_moment = desc.slrealizer.first_moment(image2)
    first_moment_x, first_moment_y = x_min+distance*first_moment[0], y_min+distance*first_moment[1]
    second_moment = desc.slrealizer.second_moment(image2)
    covariance_matrix = second_moment
    rv = scipy.stats.multivariate_normal([first_moment_x,first_moment_y], covariance_matrix, allow_singular=True) #FIX BUG     
    image = [[0]*number_of_rows for _ in range(number_of_columns)]
    image = image + rv.pdf(pos)*zeroth_moment
    print('**************zeroth moment: ', zeroth_moment)
    print('**************first moment: ', first_moment_x, first_moment_y)
    print('**************second moment: ', covariance_matrix)
    return image