#!/usr/bin/python

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from Bio import AlignIO

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

### draw a Manhattan plot of sitewise PPs and confidence windows
def manhattanPlot(df, outPath, thresholdsPP=(0, 0.8, 0.9, 1), alpha=None, beta=None,
                  xLabel="amino acid site", xSpacing=10, blkBkgd=False, width=24, height=8):

    # DataFrame column names for certain data
    # PCOC posterior probability (bar height)
    pcocPP = "PP_Max"
    # PP lower thresholds for Type I error
    typeIPP = "PP_Threshold_TypeI"
    # PP upper threshold for Type II error
    typeIIPP = "PP_Threshold_TypeII"
    # if the above columns are found then args alpha and beta will be checked for the corresponding error rates
    # if typeIIPP < typeIPP then there are insufficient data to achieve the specified error rates at that site
    # force column names for vectorized functions
    df.rename(columns={typeIPP: typeIPP, typeIIPP: typeIIPP}, inplace=True)

    # set the mode for significance calling
    # True uses sitewise sims to control error rates
    # False uses fixed PP thresholds
    simErrorControl = typeIPP in df.columns and typeIIPP in df.columns and alpha and beta

    if blkBkgd: plt.style.use('dark_background')

    # init plot
    fig, ax = plt.subplots()
    x = [i + 1 for i in range(len(df))]
    y = np.array(df[pcocPP])
    xlimits = (min(x) - 0.5, max(x) + 0.5)
    plt.xlim(xlimits)
    # set up tick interval
    xTicks = [1] + range(xSpacing, max(x), xSpacing)
    plt.xticks(xTicks, xTicks, fontsize=8, horizontalalignment='center')
    plt.yticks(fontsize=8)
    # label axes
    ax.set_xlabel(xLabel, fontsize=20)
    ax.set_ylabel("PCOC PP", fontsize=20)


    if simErrorControl:
        df = breakOutThresholds(df, alpha, beta)

        # base colors for confidence windows
        basecolorAlpha = (0.0, 0.5, 0.0) #green
        basecolorBeta = (0.5, 0.0, 0.0) #red

        # for each set of confidence thresholds
        # (len(alpha) should = len(beta))
        for i in range(len(alpha)):
            # extract pair of confidence thresholds
            ymin = np.array(df["alpha_"+str(alpha[i])])
            ymax = np.array(df["beta_"+str(beta[i])])
            # set the shade of the plot color
            # comes within 0.8 of black
            shader = np.array((1, 1, 1) * np.array([(1 - max(basecolorAlpha)) * 0.8 / len(alpha) * i] * 3))
            # make masks for sites with and without enough data

            # and plot them
            plt.bar(x, bottom=ymin, height=(1 - ymin), width=0.8, color=tuple(np.array(basecolorAlpha) + shader), alpha=0.6, zorder=1)
            plt.bar(x, bottom=0, height=ymax, width=0.8, color=tuple(np.array(basecolorBeta) + shader), alpha = 0.6, zorder=2)

        # plot PCOC PPs
        plt.bar(x, y, width=0.5, color="black", zorder=3)
        #plt.scatter(x, y, s=400, color="black", zorder=3)

    else:
        # constant threshold mode
        if not blkBkgd:
            colors = ["black", "red", "violet"]
        else:
            colors = ["white", "red", "violet"]
        # plot horizontal PP thresholds
        plt.hlines(thresholdsPP[:-1], xlimits[0], xlimits[1], colors=colors)

        # for lin
        # masks = [[pThresholds[i] >= el > pThresholds[i + 1] for el in y] for i in range(len(pThresholds) - 1)]
        # for log
        masks = np.array([[thresholdsPP[i] <= el < thresholdsPP[i + 1] for el in y] for i in range(len(thresholdsPP) - 1)])
        # plot the bars in each color
        for i, mask in enumerate(masks):
            plt.bar(x * mask, y * mask, color=colors[i], width=0.7)

    # save the plot
    fig.set_size_inches(width, height)
    plt.savefig(outPath)

    return outPath

### Draw an amino acid alignment with called columns and their trait cutoffs highlighted
### Number of red asterisks drawn on the cutoff indicates level of significance
def alignmentHighlighted(df, ali, tipTraits, outPath, highlightInconclusive = False, width=24, height=8):

    fig, ax = ali2plt(ali)

    # save the plot
    fig.set_size_inches(width, height)
    plt.savefig(outPath)

    return outPath

### Use a BioPython MSA object and embedded color scheme to plot an amino acid alignment
### Uses plt.pcolormesh()
def ali2plt(ali):

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
        'P':  (220, 150, 130)
    }

    # convert alignment object to array
    aliArray = np.array(ali)
    print aliArray
    # generate color mask
    #colorMask = np.apply(shapely.__getitem__, aliArray)
    #colorMask = np.vectorize(shapely.__getitem__)(aliArray)
    colorMask = [[None] * aliArray.shape[1]] * aliArray.shape[0]
    # doing this manually because vectorize trips on the returned tuples
    for i, row in enumerate(aliArray):
        for j, el in enumerate(row):
            colorMask[i][j] = shapely[el]
    print colorMask

    # init plot
    fig, ax = plt.subplots()
    # draw alignment
    plt.pcolormesh(aliArray)

    return fig, ax


### Helper function to unpack lists of confidence thresholds into their own columns in the dataframe
def breakOutThresholds(df, alpha, beta):
    # split out the significance threshold lists
    df[["alpha_" + str(thres) for thres in alpha]] = pd.DataFrame(df.PP_Threshold_TypeI.tolist(), index=df.index)
    df[["beta_" + str(thres) for thres in beta]] = pd.DataFrame(df.PP_Threshold_TypeII.tolist(), index=df.index)

    return df