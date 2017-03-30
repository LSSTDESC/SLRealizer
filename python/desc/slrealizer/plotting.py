import plotting
import pandas as pd
import matplotlib.pyplot as plt
import pylab
import matplotlib
import math

def plot_lens_random_date(self, lensID = 7176527, convolve=False):
    import random
    randomIndex = random.randint(0, 200)
    if(self.observation[randomIndex][1]!='y'):
        currLens = self.lens.get_lens(lensID)
        drawPlot(self.observation[randomIndex], currLens, convolve)
    else:
        # recursive call
        plot_lens_random_date(self, lensID, convolve)

def determineAlpha(mag_ratio):
    if(mag_ratio>1):
        quasar_alpha = 1/mag_ratio
        lens_alpha = 1
    else:
        quasar_alpha = 1
        lens_alpha = mag_ratio
    return quasar_alpha, lens_alpha

def determineScale(lensX, sourceX, lensY, sourceY):
    if (abs(lensX[0]+sourceX))>(abs(sourceX)):
        plotX = abs(lensX[0]+sourceX)
    else:
        plotX = abs(sourceX)
    if (abs(lensY[0]+sourceY))>(abs(sourceY)):
        plotY = abs(lensY[0]+sourceY)
    else:
        plotY = abs(sourceY)
    return plotX, plotY

def plotFigureOnMatplotlib(currObs, convolve, quasar_alpha, lens_alpha, sourceX, sourceY, lensX, lensY, plotX, plotY):
    #fig = plt.figure(1)
    scale_factor = 2
    if(currObs[1]=='g'):
        circleColor = 'g'
    else:
        circleColor = 'r'
    plt.axis('scaled')
    plt.ylim(-3, 3)
    plt.xlim(-3, 3)
    plt.xlabel('xPosition')
    plt.ylabel('yPosition')
    plt.title('Observation with ' + 'filter ' + currObs[1] + ' on MJD ' + str(currObs[0]))
    for i in range(4):
        source = plt.Circle((sourceX, sourceY), radius=currObs[2]/scale_factor, alpha=quasar_alpha, fc=circleColor, linewidth=0)
        if(convolve):
            # For now, we just let lensFWHM to be zero.
            #lensFWHM = currLens['REFF']
            lensFWHM = 0.0
            lens = plt.Circle((lensX[i]+sourceX, lensY[i]+sourceY), radius=convolve(currObs[2], lensFWHM)/scale_factor, alpha=lens_alpha, fc=circleColor, linewidth=0)
        else:
            lens = plt.Circle((lensX[i]+sourceX, lensY[i]+sourceY), radius=currObs[2]/scale_factor, alpha=lens_alpha, fc=circleColor, linewidth=0, edgecolor='b', linestyle='dashed')
        plt.gca().add_patch(lens)
        plt.gca().add_patch(source)

        #seeing = plt.Circle(((-plotY-2*currObs[2])/scale_factor, (plotX-2*currObs[2])/scale_factor), radius=currObs[2]/scale_factor, alpha=0.1, fc='k')
    seeing = plt.Circle((-2.5, -2.5), radius=currObs[2]/scale_factor, alpha=0.1, fc='k')
    plt.legend((source, seeing, lens),('source', 'seeing', 'lens'))
    plt.gca().add_patch(seeing)


def drawPlot(currObs, currLens, convolve = False):
    #obsHist has MJD Filter FWHM 5sigmag
    filterQuasar = currObs[1] + '_SDSS_quasar'
    filterLens = currObs[1] + '_SDSS_lens'
    lens_mag = currLens[filterLens]
    quasar_mag = currLens[filterQuasar]
    mag_ratio = math.pow(2.5, -lens_mag+quasar_mag)
    quasar_alpha, lens_alpha = determineAlpha(mag_ratio)
    sourceX = currLens['XSRC'][0]
    sourceY = currLens['YSRC'][0]
    lensX = currLens['XIMG'][0]
    lensY = currLens['YIMG'][0]
    plotX, plotY = determineScale(lensX, sourceX, lensY, sourceY)
    plotFigureOnMatplotlib(currObs, convolve, quasar_alpha, lens_alpha, sourceX, sourceY, lensX, lensY, plotX, plotY)

def convolve(obsFWHM, initialFWHM = 0.0):
    seeing = obsFWHM
    seeingsigma = seeingToSig(seeing)
    initSigma = seeingToSig(initialFWHM)
    convolveSigma = seeingsigma + initSigma
    return convolveSigma
    
def seeingToSig(seeing):
    return seeing / np.sqrt(8 * np.log(2))

def sigmaTofwhm(sigma):
    return sigma * np.sqrt(8 * np.log(2))
