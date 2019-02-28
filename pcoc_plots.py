#!/usr/bin/python

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors

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

### Draw heatmap with a wide pd.DataFrame() with indices to be used as labels, an output path and optionally:
### an alignment to draw seq logo over the heatmap
### a colormap to draw a colorimetric legend along the y-axis.
def heatMapDF(df, outfile, rainbow=False):

    cutoffs = df.index.values.tolist()
    sites = df.columns.values.tolist()

    pcocPP = df.values

    # for dark bkgd
    plt.style.use('dark_background')
    fig, ax = plt.subplots()

    if rainbow:
        # add color column to heatmap
        sites = ["color"] + sites
        # make a mask for all but the first column
        mask = np.array(pcocPP.shape[0] * [[False] + ([True] * pcocPP.shape[1])], dtype=bool)
        # add cutoffs to PP array as first column
        pcocPP = np.hstack((np.array([cutoffs]).T, pcocPP)) # populate that column with cutoff vals

        colorscale = np.ma.masked_where(mask, pcocPP)
        heatmap = np.ma.masked_where(np.invert(mask), pcocPP)

        # get the custom colormap and normalization
        rbMap, rbNorm = realRainbow()
        ax.pcolormesh(colorscale, cmap=rbMap, norm=rbNorm)
        ppMesh = ax.pcolormesh(heatmap, vmin=0., vmax=1., cmap='Greys')
    else:
        ppMesh = ax.pcolormesh(pcocPP, vmin=0., vmax=1., cmap='Greys')

    # for PCOC gene results
    # Figure labels
    ax.set_xlabel("site", fontsize=20)
    plt.xticks(np.arange(len(sites)) + 0.5, ["color"] + [int(float(site)) for site in sites[1:]], rotation=90, size=7)
    ax.xaxis.set_label_position('top')
    ax.xaxis.set_ticks_position('top')
    ax.set_ylabel("detection cutoff for emission wavelength (nm)", fontsize=20)
    plt.yticks(np.arange(len(cutoffs)) + 0.5, [int(cutoff) for cutoff in cutoffs], size=20)
    ax.invert_yaxis()
    ax.patch.set_facecolor('white')
    # for spine in ax.spines: ax.spines[spine].set_color('white')
    cbar = fig.colorbar(ppMesh)
    cbar.ax.text(0.55, 0.025, 'PCOC PP', rotation=90, ha='center', va='bottom', transform=cbar.ax.transAxes,
                 fontsize=8, color='black')

    fig.tight_layout()
    # plt.show()
    fig.set_size_inches(24, 8)
    #plt.savefig(outfile, transparent=True)
    plt.savefig(outfile)

    '''
    # for PCOC calling-accuracy heatmaps
    # Figure labels
    ax.set_xlabel("detection cutoff (nm)", fontsize=20)
    plt.xticks(np.arange(len(sites)) + 0.5, ["color"] + [int(float(site)) for site in sites[1:]], rotation=90, size=20)
    ax.xaxis.set_label_position('bottom')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_ylabel("simulation cutoff (nm)", fontsize=20)
    plt.yticks(np.arange(len(cutoffs)) + 0.5, [int(cutoff) for cutoff in cutoffs], size=20)
    ax.invert_yaxis()
    ax.patch.set_facecolor('white')
    #for spine in ax.spines: ax.spines[spine].set_color('white')
    cbar = fig.colorbar(ppMesh)
    cbar.ax.text(0.55, 0.025, 'calling rate', rotation=90, ha='center', va='bottom', transform=cbar.ax.transAxes,
                 fontsize=20, color='black')

    fig.tight_layout()
    # plt.show()
    fig.set_size_inches(12, 10)
    plt.savefig(outfile, transparent=True)
    '''

# a command line interface to make a heatmap from tsv files. If a directory full of files is passed,
# all files in that dir will be averaged elementwise
def main(argv, wayout):
    if not len(argv):
        argv.append('-h')
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    parser.add_argument("-i", "--input", type=str, help="source of tsv file(s)", required=True)
    parser.add_argument("-o", "--output", type=str, help="name of PDF output", required=True)

    global args
    args = parser.parse_args(argv)

    # HFS+ - safe listdir() function
    def listdirSafe(path):
        # return sorted([path + '/' + f for f in os.listdir(path) if f[:2] != '._'])
        return sorted([f for f in os.listdir(path) if f[:1] != '.'])

    # take a list of DataFrames and average them elementwise
    # use indices from the first df in the list
    # requires numpy
    def dfMean(dfList):
        rows = dfList[0].index
        cols = dfList[0].columns
        matrixList = [df.values for df in dfList]

        # sum the matrices elementwise
        outMatrix = np.zeros_like(matrixList[0])
        for matrix in matrixList:
            outMatrix = np.add(outMatrix, matrix)

        # divide to get the mean
        outMatrix = outMatrix / len(matrixList)

        # convert back to df and name the axes
        outDf = pd.DataFrame(outMatrix)
        outDf.index = rows
        outDf.columns = cols

        return outDf

    # make heatmaps
    if os.path.isdir(args.input):
        # load all the dfs in the directory
        infiles = [args.input + '/' + infile for infile in listdirSafe(args.input)]
        dfs = [pd.read_table(infile, index_col=0) for infile in infiles]
        # get elementwise mean df
        meanDf = dfMean(dfs)
        # draw heatmap
        heatMapDF(meanDf, args.output, rainbow=True)
        # print the averaged df
        print >> wayout, meanDf.to_csv(sep='\t')
    else:
        heatMapDF(pd.read_table(args.input), args.output, rainbow=True)


### draw a Manhattan plot of the summed PPs
def manhattanPlot(df, outPath, keyID=None):
    # at present this is meant to be run once
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    x = [i + 1 for i in range(len(df))]
    y = np.array(df)
    # y = np.array([-math.log10(p) for p in y]) # log-transform
    xlimits = (min(x) - 0.5, max(x) + 0.5)
    plt.xlim(xlimits)
    plt.xticks(x, x, rotation=90, fontsize=7)
    plt.yticks(fontsize=20)

    if keyID:
        ax.set_xlabel("position in " + keyID, fontsize=20)
    else:
        ax.set_xlabel("alignment column")
    ax.set_ylabel("total PCOC PP", fontsize=20)

    # pThresholds = (0, 0.8, 0.9, 0.95, 1)
    # colors = ["black", "blue", "red", "violet"]
    pThresholds = (0, 0.8, 1)
    colors = ["white", (0, 1, 0.05)]  # for black bkgd
    # colors = ["black", "blue"]  # for white bkgd
    plt.hlines(pThresholds[:-1], xlimits[0], xlimits[1], colors=colors)

    # for lin
    # masks = [[pThresholds[i] >= el > pThresholds[i + 1] for el in y] for i in range(len(pThresholds) - 1)]
    # for log
    masks = np.array([[pThresholds[i] <= el < pThresholds[i + 1] for el in y] for i in range(len(pThresholds) - 1)])

    for i, mask in enumerate(masks):
        plt.bar(x * mask, y * mask, color=colors[i], width=1)

    fig.set_size_inches(24, 8)
    # plt.show()
    plt.savefig(outPath)

    # no return

if __name__ == "__main__":
    # when called from command line, print results to stdout
    main(sys.argv[1:], sys.stdout)