#!/usr/bin/python

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
from Bio import AlignIO

### draw PCOC.ontinuous master figure with the specified elements
### 0 plots everything
### 1 plots just the Manhattan
### 2 plots just the highlighted alignment
def masterFigure(df, ali, tipTraits, elements=0, alpha=None, beta=None, thresPP=[0.8, 0.9, 0.95], blkBkgd=False,
                 width=24, height=6, xLabel="amino acid site", fontsize=None, outPath=None):

    if blkBkgd: plt.style.use('dark_background')

    if elements == 0:

        if not fontsize:
            fontsize = min(2 * height, 75 * width / len(df))

        # set up axes w/equal size
        #fig, ax = plt.subplots(2, 1)
        # set up axes w/adjustable height ratio
        fig = plt.figure()
        spec = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[1,2])
        ax = [fig.add_subplot(spec[row, 0]) for row in range(2)]

        # plot Manhattan
        plt.subplot(ax[0])
        manhattanPlot(df, alpha=alpha, beta=beta, thresPP=thresPP, fontsize=fontsize)
        # take away the xlabel
        plt.xlabel('')
        # plot alignment
        plt.subplot(ax[1])
        alignmentHighlighted(df, ali, tipTraits, xLabel=xLabel, fontsize=fontsize)


    elif elements == 1:
        if not fontsize:
            fontsize = min(4 * height, 75 * width / len(df))

        fig, ax = plt.subplots()
        # plot Manhattan
        manhattanPlot(df, alpha=alpha, beta=beta, thresPP=thresPP, xLabel=xLabel, fontsize=fontsize)

    elif elements == 2:
        if not fontsize:
            fontsize = min(4 * height, 75 * width / len(df))

        fig, ax = plt.subplots()
        alignmentHighlighted(df, ali, tipTraits, xLabel=xLabel, fontsize=fontsize)

    # save the figure
    plt.tight_layout()
    fig.set_size_inches(width, height)
    plt.savefig(outPath)



### draw a Manhattan plot of sitewise PPs and confidence windows
def manhattanPlot(df, thresPP=[0.8, 0.9, 0.95], alpha=None, beta=None,
                  xLabel="amino acid site", xSpacing=10, blkBkgd=False, width=20, height=8, fontsize=None, outPath=None):

    # DataFrame column names for certain data
    # PCOC posterior probability (bar height)
    pcocPP = "PP_Max"
    # PP lower thresholds for Type I error
    typeIPP = "PP_Threshold_TypeI"
    # PP upper thresholds for Type II error
    typeIIPP = "PP_Threshold_TypeII"
    # categorical significance level
    sigLevel = "SigLevel"
    # if the above columns are found then args alpha and beta will be checked for the corresponding error rates
    # if typeIIPP < typeIPP then there are insufficient data to achieve the specified error rates at that site
    # force column names for vectorized functions
    df.rename(columns={typeIPP: typeIPP, typeIIPP: typeIIPP}, inplace=True)
    # SORTING IS IMPORTANT DOWNSTREAM
    alpha = sorted(alpha, reverse=True)
    beta = sorted(beta)
    thresPP = [0] + thresPP

    # if font size not specified by caller, attempt to autoscale
    if not fontsize:
        fontsize = min(2.5 * height, 75 * width / len(df))

    # set the mode for significance calling
    # True uses sitewise sims to control error rates
    # False uses fixed PP thresholds
    simErrorControl = typeIPP in df.columns and typeIIPP in df.columns and alpha and beta

    if blkBkgd: plt.style.use('dark_background')

    # init plot
    x = [i + 1 for i in range(len(df))]
    y = np.array(df[pcocPP])
    xlimits = (min(x) - 0.5, max(x) + 0.5)
    plt.xlim(xlimits)
    # set up tick interval
    xTicks = [1] + range(xSpacing, max(x), xSpacing)
    plt.xticks(xTicks, xTicks, fontsize=fontsize*2, horizontalalignment='center')
    plt.yticks(fontsize=fontsize*2)
    # label axes
    plt.xlabel(xLabel, fontsize=fontsize*4)
    plt.ylabel("PCOC PP", fontsize=fontsize*4)


    if simErrorControl:
        df = breakOutThresholds(df, alpha, beta)

        # base colors for confidence windows
        basecolorAlpha = (0.0, 0.5, 0.0) #green
        basecolorBeta = (0.5, 0.0, 0.0) #red
        # alpha for the lowest level of significance
        baseSolidness = 0.6

        # for each set of confidence thresholds
        for i in range(len(alpha)):
            # extract pair of confidence thresholds
            ymin = np.array(df["alpha_"+str(alpha[i])])
            #ymin = np.array(df["alpha_" + str(i)])

            # and plot them
            solidness = baseSolidness + i/len(alpha) * (1 - baseSolidness)
            plt.bar(x, bottom=ymin, height=(1 - ymin), width=0.8, color=basecolorAlpha, alpha=solidness)

        for i in range(len(beta)):
            # extract pair of confidence thresholds
            ymax = np.array(df["beta_"+str(beta[i])])
            #ymax = np.array(df["beta_" + str(i)])

            # and plot them
            solidness = baseSolidness + i/len(beta) * (1 - baseSolidness)
            plt.bar(x, bottom=0, height=ymax, width=0.8, color=basecolorBeta, alpha=solidness)

        # plot PCOC PPs
        ax = plt.gca()
        #xUnit = 12
        # width of a single unit in pixels
        xUnit = (ax.transData.transform((1, 0)) - ax.transData.transform((0, 0)))[0]
        #print >> sys.stderr, "1 x unit={}".format(xUnit)

        #pointArea = 200.0/df.shape[0]
        #if pointArea < 1:
        #    pointArea = 1
        #else:
        #    pointArea = 2.25*float(round(pointArea))
        plt.scatter(x, y, s=8*xUnit**2, marker="D", color="black", zorder=3)  # s = area of 400 pt^2

    else:
        # constant threshold mode
        if not blkBkgd:
            colors = ["black", "blue", "red", "violet"]
        else:
            colors = ["white", "blue", "red", "violet"]

        # plot horizontal PP thresholds
        plt.hlines(thresPP[1:], xlimits[0], xlimits[1], colors=colors[1:])

        # for lin
        # masks = [[pThresholds[i] >= el > pThresholds[i + 1] for el in y] for i in range(len(pThresholds) - 1)]
        # for log

        # old mask compares PP with thresholds
        #masks = np.array([[thresPP[i] <= el < thresPP[i + 1] for el in y] for i in range(len(thresPP) - 1)])

        # new mask just checks the passed sig level
        masks = np.array([[el == i for el in df[sigLevel]] for i in range(len(thresPP))])

        # plot the bars in each color
        for i, mask in enumerate(masks):
            plt.bar(x * mask, y * mask, color=colors[i], width=0.7)

    # save the plot
    if outPath:
        plt.figure(figsize=(width, height))
        plt.savefig(outPath)

### Draw an amino acid alignment with called columns and their trait cutoffs highlighted
### Number of red asterisks drawn on the cutoff indicates level of significance
def alignmentHighlighted(df, ali, tipTraits, xLabel="amino acid site", blkBkgd=False, width=20, height=8, fontsize=None, outPath=None):

    if blkBkgd: plt.style.use('dark_background')

    # DataFrame column names for certain data
    # PCOC posterior probability (bar height)
    pcocPP = "PP_Max"
    # categorical significance level
    sigLevel = "SigLevel"
    # trait cutoff
    cutoff = "Cutoff"
    '''
    # deprecate these - 20190313
    # PP lower thresholds for Type I error
    typeIPP = "PP_Threshold_TypeI"
    # PP upper threshold for Type II error
    typeIIPP = "PP_Threshold_TypeII"
    '''
    # if font size not specified by caller, attempt to autoscale
    if not fontsize:
        fontsize = min(2.5 * height, 75 * width / len(df))

    # reorder alignment according to trait value
    ali.sort(key = lambda record: tipTraits[record.id])

    # cull the tipTraits dict so it only contains taxa found in ali
    alignTaxa = [record.id for record in ali]
    tipTraits = {key:value for key, value in tipTraits.items() if key in alignTaxa}

    # plot basic MSA
    ali2plt(ali, fontsize=fontsize, rowSpacing=fontsize/10)
    plt.xlabel(xLabel, fontsize=fontsize*4)

    # push this back to the parent script and have it fed to the plotting routine
    # add a significance level column to the dataFrame: 0 (not sig) -> index of the most stringent alpha cutoff
    #df[sigLevel] = np.vectorize(rank)(df[pcocPP], df[typeIPP])

    # set up a colormap that varies alpha
    cm = colors.LinearSegmentedColormap.from_list(name = 'cm', colors = [(1,1,1,0.8),(1,1,1,0)])

    # plot a colormesh on top of the alignment axes with alpha mapped to sigLevel
    greyOut = np.array([df[sigLevel].tolist()] * len(ali))
    plt.pcolormesh(greyOut, cmap=cm, norm=colors.Normalize(vmin=0, vmax=max(df[sigLevel]), clip=False), zorder = 3)
    #fig.colorbar(im, ax=ax)

    # draw a patch with asterisks on every significant cutoff
    ax = plt.gca()
    xUnit = (ax.transData.transform((1, 0)) - ax.transData.transform((0, 0)))[0]
    x = [s + 0.5 for s in range(df.shape[0]) if df[sigLevel].tolist()[s] > 0]
    sigCuts = [df[cutoff].tolist()[j] for j in range(df.shape[0]) if df[sigLevel].tolist()[j] > 0]
    sortCuts = sorted(tipTraits.values())
    y = [rank(cut, sortCuts) for cut in sigCuts]
    ## custom marker shape
    rectMarker = [(-2, -0.5), (-2, 0.5), (2, 0.5), (2, -0.5)]
    plt.scatter(x, y, s=3 * (xUnit ** 2), marker=rectMarker, color="black", zorder=4)
    plt.scatter(x, y, s=(xUnit ** 2) / 16, marker="*", color="red", zorder=5)

    # add triangle markers at top of alignment
    # x-coords are same as above
    #x = [s + 0.5 for s in range(df.shape[0]) if df[sigLevel].tolist()[s] > 0]
    ### push ylim up a little
    plt.ylim((plt.ylim()[0], plt.ylim()[1] - 2))
    ### to accommodate this
    y = [plt.ylim()[1] + 0.5] * len(x)
    plt.scatter(x, y, s=32 * xUnit ** 2, marker="v", color="red", zorder=6)  # s = area of 400 pt^2
    ax.spines['top'].set_visible(False)
    # fix x axis
    plt.xlim((0, df.shape[0]))


    # save the plot
    if outPath:
        plt.figure(figsize=(width, height))
        plt.savefig(outPath)

    #return fig, ax

### Use a BioPython MSA object and embedded color scheme to plot an amino acid alignment
### Uses plt.pcolormesh()
def ali2plt(ali, rowSpacing=None, xSpacing=10, yLabel=None, blkBkgd=False, fontsize=None):

    # color schemes
    shapely = {
        'D':  (230, 10, 10),
        'E':  (230, 10, 10),
        'C':  (230, 230, 0),
        'M':  (230, 230, 0),
        'K':  (20, 90, 255),
        'R':  (20, 90, 255),
        'S':  (250, 150, 0),
        'T':  (250, 150, 0),
        'F':  (50, 50, 170),
        'Y':  (50, 50, 170),
        'N':  (0, 220, 220),
        'Q':  (0, 220, 220),
        'G':  (235, 235, 235),
        'L':  (15, 130, 15),
        'V':  (15, 130, 15),
        'I':  (15, 130, 15),
        'A':  (200, 200, 200),
        'W':  (180, 90, 180),
        'H':  (130, 130, 210),
        'P':  (220, 150, 130),
        '-':  (255, 255, 255)
    }
    # scale it down to [0,1] RGB
    for key, value in shapely.iteritems():
        shapely[key] = tuple([x/255.0 for x in value])

    # turn dict into a colormap
    # first get the letters, integer indices and colors into parallel lists
    letters = [None] * 21
    ints = [None] * 21
    chroma = [None] * 21
    ## CHANGE COLOR SCHEME HERE ##
    for i, key in enumerate(sorted(shapely.keys())):
        letters[i] = key
        ints[i] = i
        chroma[i] = shapely[key]
    #cm = colors.LinearSegmentedColormap.from_list("protein", chroma, N=20)
    # This colormap is not behaving accurately!
    cm = colors.ListedColormap(chroma)

    numSeqs = len(ali)
    numSites = ali.get_alignment_length()
    seqNames = [record.id for record in ali]

    # if font size not specified by caller, attempt to autoscale
    if not fontsize:
        fontsize = 180 / ali.get_alignment_length()

    # convert strings to ints for colormapping
    aliArray = np.array(ali)
    intdices = {key: value for key, value in zip(letters, ints)}
    cellData = np.vectorize(intdices.__getitem__)(aliArray)

    # x-axis formatting
    xTicks = [1] + range(xSpacing, numSites-1, xSpacing)
    plt.xticks([x - 0.5 for x in xTicks], xTicks, fontsize=fontsize*2, horizontalalignment='center')

    # y-axis formatting
    # This is needed to make the first row of the array correspond to 1st seq in MSA
    plt.gca().invert_yaxis()
    plt.yticks([y + 0.5 for y in range(numSeqs)], seqNames, fontsize=fontsize*2, verticalalignment='center')
    #ax.axes.get_yaxis().set_visible(False)

    #TODO: use pcolormesh(X, Y) to put breathing room between the aligned sequences

    # draw alignment
    plt.pcolormesh(cellData, cmap=cm, norm=colors.Normalize(vmin=-0.5, vmax=20.5, clip=False), zorder = 1)
    # to check proper normalization of colorscale
    # fig.colorbar(im, ax=ax)

    # ylabel
    if not yLabel:
        plt.ylabel('species', fontsize=fontsize*4)

    # put letters on the alignment
    for x in range(numSites):
        for y in range(numSeqs):
            plt.annotate(str(aliArray[y][x]), (x + 0.5, y + 0.5), fontsize=fontsize, horizontalalignment='center', verticalalignment='center', zorder=2)

    # give the alignment rows breathing room
    if rowSpacing:
        if blkBkgd: chrom = "black"
        else: chrom = "white"
        plt.hlines([0] + range(numSeqs), 0, numSites, colors=chrom, linewidth=rowSpacing)

## HELPER FUNCTIONS ##

### Helper function to unpack lists of confidence thresholds into their own columns in the dataframe
def breakOutThresholds(df, alpha, beta):

    if isinstance(alpha, float): alpha = [alpha]
    if isinstance(alpha, float): beta = [beta]
    # split out the significance threshold lists
    for i, thres in enumerate(alpha):
        df["alpha_" + str(thres)] = [j[i] if isinstance(j, list) and len(j) else np.nan for j in df["PP_Threshold_TypeI"].tolist()]
    for i, thres in enumerate(beta):
        df["beta_" + str(thres)] = [j[i] if isinstance(j, list) and len(j) else np.nan for j in df["PP_Threshold_TypeII"].tolist()]

    return df

### Rank val in iterable seri
def rank(val, seri):
    if isinstance(seri, list):
        rank = 0
        for i, thres in enumerate(sorted(seri)):
            if val >= thres:
                rank += 1
        return rank
    else:
        return 0

### Create a discrete colormap for monochromatic light based on wavelength
### FPs are not monochromatic, so my work is not done here!
def realRainbow():

    import wl2rgb

    # colormap should go from 380-750 nm
    wls = range(380, 750)
    chroma = [wl2rgb.wavelength_to_rgb(wl, gamma=0.8, scaleMax=1) for wl in wls]
    cmap = colors.ListedColormap(chroma)
    boundaries = [min(wls)-1] + wls # fenceposted!
    norm = colors.BoundaryNorm(boundaries, cmap.N, clip=True)
    return cmap, norm