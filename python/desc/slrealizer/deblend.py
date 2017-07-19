from scipy.stats import multivariate_normal
import plotting
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, MaxNLocator
from numpy import linspace
import matplotlib
import scipy
import scipy.optimize as opt
import math

"""
    
    deblend.py deblends the light sources for OM10. This assumes null deblender where all the sources are observed as a single object. All the sources are assumed to have Gaussian PSFs.
    """

def deblend(currObs, currLens):
    print 'Deblending starts.....'
    print 'The system looks like this....'
    #plotting.draw_model(currObs, currLens)
    blend_all_objects(currObs, currLens)

def blend_all_objects(currObs, currLens):
    print 'printing currLens'
    #obsHist has MJD Filter FWHM 5sigmag
    filterQuasar = currObs[1] + '_SDSS_quasar'
    filterLens = currObs[1] + '_SDSS_lens'
    lens_mag = currLens[filterLens]
    quasar_mag = currLens[filterQuasar]
    mag_ratio = math.pow(2.5, -lens_mag+quasar_mag)
    sourceX = currLens['XSRC'][0]
    sourceY = currLens['YSRC'][0]
    x_min = -3
    x_max = 3
    y_min = -3
    y_max = 3
    distance = 0.01
    x, y = np.mgrid[x_min:x_max:distance, y_min:y_max:distance]
    pos = np.dstack((x, y))
    PSF_HWHM = currObs[2]
    number_of_rows = int((x_max - x_min)/distance)
    number_of_columns = int((y_max - y_min)/distance)
    rv = scipy.stats.multivariate_normal([sourceX,sourceY], [[PSF_HWHM*PSF_HWHM, 0], [0, PSF_HWHM*PSF_HWHM]])
    image = [[0]*number_of_rows for _ in range(number_of_columns)]
    image = image + rv.pdf(pos)
    print 'PSF_HWHM: ',
    print PSF_HWHM
    for i in xrange(4):
        print currLens['XIMG'][0][i]
        print currLens['YIMG'][0][i]
        rv = scipy.stats.multivariate_normal([currLens['XIMG'][0][i],currLens['YIMG'][0][i]], [[PSF_HWHM*PSF_HWHM, 0], [0, PSF_HWHM*PSF_HWHM]])
        image = image + rv.pdf(pos)
    cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list('my_colormap', ['black', 'green', 'yellow'], 256)
    img2 = plt.imshow(image, interpolation='nearest', cmap = cmap2, origin='lower', extent=[-3, 3, -3, 3], aspect = "auto")
    plt.colorbar(img2,cmap=cmap2)
    plt.show()

    import photutils
    segm = photutils.detect_sources(image, 0.1, 3)
    # plot the image and the segmentation image
    plt.imshow(segm.data, origin='lower', interpolation='nearest')

    from photutils import DAOStarFinder

    daofind = DAOStarFinder(fwhm=1.5, threshold=0.00001)
    sources = daofind(image)
    from astropy.visualization import SqrtStretch
    from astropy.visualization.mpl_normalize import ImageNormalize
    from photutils import CircularAperture
    positions = (sources['xcentroid'], sources['ycentroid'])
    apertures = CircularAperture(positions, r=4.)
    norm = ImageNormalize(stretch=SqrtStretch())
    plt.imshow(image, cmap='Greys', origin='lower', norm=norm)
    apertures.plot(color='red', lw=3.0, alpha=1.0)
    print sources



#https://stackoverflow.com/questions/21566379/fitting-a-2d-gaussian-function-using-scipy-optimize-curve-fit-valueerror-and-m
