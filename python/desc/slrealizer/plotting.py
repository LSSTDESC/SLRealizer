fig = pylab.figure()
scale_factor = 10
seeing = plt.Circle((0.8, 0.8), radius=obsHist[1][2]/scale_factor, alpha=0.5, fc='grey')
plt.gca().add_patch(seeing)