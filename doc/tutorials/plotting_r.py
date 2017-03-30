import pylab
import matplotlib
import math
import matplotlib.pyplot as plt

def plot_lens(lens, obsHist, convolve=False):
    #obsHist has MJD Filter FWHM 5sigmag
    filterQuasar = obsHist[1] + '_SDSS_quasar'
    filterLens = obsHist[1] + '_SDSS_lens'
    lens_mag = lens[filterLens]
    quasar_mag = lens[filterQuasar]
    mag_ratio = math.pow(2.5, -lens_mag+quasar_mag)
    if(mag_ratio>1):
        quasar_alpha = 1/mag_ratio
        lens_alpha = 1
    else:
        quasar_alpha = 1
        lens_alpha = mag_ratio
    scale_factor = 2
    sourceX = lens['XSRC'][0]
    sourceY = lens['YSRC'][0]
    lensX = lens['XIMG'][0]
    lensY = lens['YIMG'][0]
    if (abs(lensX[0]+sourceX))>(abs(sourceX)):
        plotX = abs(lensX[0]+sourceX)
    else:
        plotX = abs(sourceX)
    if (abs(lensY[0]+sourceY))>(abs(sourceY)):
        plotY = abs(lensY[0]+sourceY)
    else:
        plotY = abs(sourceY)
    fig = plt.figure(figsize=(4, 16))

    if(True):
        for i in range(4):
            if(obsHist[1]=='g'):
                circleColor = 'g'
            else:
                circleColor = 'r'
            plotNumStr = "4"+"1"+str(i+1)
            plotNum = int(plotNumStr)
            sub = fig.add_subplot(plotNum)
            sub.set_ylim(-plotY-2*obsHist[2]/scale_factor, plotY+2*obsHist[2]/scale_factor)
            sub.set_xlim(-plotX-2*obsHist[2]/scale_factor, plotX+2*obsHist[2]/scale_factor)
            source = plt.Circle((sourceX, sourceY), radius=obsHist[2]/scale_factor, alpha=quasar_alpha, fc=circleColor, linewidth=0)
            lens = plt.Circle((lensX[i]+sourceX, lensY[i]+sourceY), radius=obsHist[2]/scale_factor, alpha=lens_alpha, fc=circleColor, linewidth=0)
            fig.gca().add_patch(lens)
            fig.gca().add_patch(source)
            sub.set_title('Observation ' + str(i+1) + ' with ' + 'filter ' + obsHist[1])
            seeing = plt.Circle(((-plotY-2*obsHist[2])/scale_factor, (plotX-2*obsHist[2])/scale_factor), radius=obsHist[2]/scale_factor, alpha=0.1, fc='k')
            sub.legend((source, seeing, lens),('source', 'seeing', 'lens'))
            plt.gca().add_patch(seeing)

