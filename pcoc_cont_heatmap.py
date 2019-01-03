#!/usr/bin/python

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

## Contains a routine for plotting site x cutoff heatmaps for output from molecular convergence analyses.
## Based on example from https://matplotlib.org/gallery/images_contours_and_fields/image_annotated_heatmap.html

def heatmap(data, row_labels, col_labels, ax=None, cax=None, cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Arguments:
        data       : A 2D numpy array of shape (N,M)
        row_labels : A list or array of length N with the labels
                     for the rows
        col_labels : A list or array of length M with the labels
                     for the columns
    Optional arguments:
        ax         : A matplotlib.axes.Axes instance to which the heatmap
                     is plotted. If not provided, use current axes or
                     create a new one.
        cbar_kw    : A dictionary with arguments to
                     :meth:`matplotlib.Figure.colorbar`.
        cbarlabel  : The label for the colorbar
    All other arguments are directly passed on to the imshow call.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    if cax is not None:
        cbar = ax.figure.colorbar(im, fraction=0.005, pad=0.01, **cbar_kw)
        cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=90, ha="left", va='center',
             rotation_mode="anchor")
    #plt.setp(ax.get_xticklabels()) #TEST

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    #ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    if cax is not None:
        return im, cbar
    else:
        return im


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=["black", "white"],
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Arguments:
        im         : The AxesImage to be labeled.
    Optional arguments:
        data       : Data used to annotate. If None, the image's data is used.
        valfmt     : The format of the annotations inside the heatmap.
                     This should either use the string format method, e.g.
                     "$ {x:.2f}", or be a :class:`matplotlib.ticker.Formatter`.
        textcolors : A list or array of two color specifications. The first is
                     used for values below a threshold, the second for those
                     above.
        threshold  : Value in data units according to which the colors from
                     textcolors are applied. If None (the default) uses the
                     middle of the colormap as separation.

    Further arguments are passed on to the created text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[im.norm(data[i, j]) > threshold])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts


def realRainbow():
    # colormap should go from 380-750 nm

    import wl2rgb
    from matplotlib.colors import LinearSegmentedColormap

    wls = range(380, 750)
    colors = [wl2rgb.wavelength_to_rgb(wl, gamma=0.8, scaleMax=1) for wl in wls]
    return LinearSegmentedColormap.from_list('int_rainbow', colors, N=len(wls))


# Expects a wide pd.DataFrame() with indices to be used as labels, an output path and optionally:
# an alignment to draw seq logo over the heatmap
# a colormap to draw a colorimetric legend along the y-axis.
def heatMapDF(df, outfile, rainbow=False):

    cutoffs = df.index.values.tolist()
    sites = df.columns.values.tolist()

    pcocPP = df.values

    #plt.style.use('dark_background')

    #heatmap(pcocPP, cutoffs, sites, ax=axs[1], cax=axs[1], cmap="gray", cbarlabel="PCOC PP")

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

        ax.pcolormesh(colorscale, cmap=realRainbow())
        ppMesh = ax.pcolormesh(heatmap, cmap='Greys')
    else:
        ppMesh = ax.pcolormesh(pcocPP, cmap='Greys')

    # Figure labels
    ax.set_xlabel("position in avicGFP", fontsize=20)
    plt.xticks(np.arange(len(sites)) + 0.5, sites, rotation=90)
    ax.xaxis.set_label_position('top')
    ax.xaxis.set_ticks_position('top')
    ax.set_ylabel("emission WL cutoff (nm)", fontsize=20)
    plt.yticks(np.arange(len(cutoffs)) + 0.5, [int(cutoff) for cutoff in cutoffs])
    ax.invert_yaxis()
    cbar = fig.colorbar(ppMesh)
    cbar.ax.text(0.55, 0.025, 'PCOC PP', rotation=90, ha='center', va='bottom', transform=cbar.ax.transAxes, fontsize=20)

    fig.tight_layout()
    #plt.show()
    fig.set_size_inches(48, 12)
    #fig.set_size_inches(70, 10)  # 70" wide @340 col algt. #TEST it's a magic number!
    plt.savefig(outfile)