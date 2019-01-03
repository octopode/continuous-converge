#!/usr/bin/env python

"""tpr-heatmapper.py
v1.0 (c)2018 Jacob Winnikoff
jwinnikoff@mbari.org

Parses file structures of converge-detect results to plot true positive, false negative rates

Input: directories containing ancestral (non-convergent) and convergent simulations
Output: tab-delimited tables of TPR, FNR for different levels of simulation and detection cutoff, heatmaps if desired

Usage:

tpr-heatmapper.py -a PCOC-anc -c PCOC-con -t 0.8 -h PCOC_sim_heatmaps.pdf > PCOC_sim_TPR-FNR.tsv
"""

import sys
import os
import argparse
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pcoc_cont_heatmap as heatmapper

# HFS+ - safe listdir() function
def listdirSafe(path):
    #return sorted([path + '/' + f for f in os.listdir(path) if f[:2] != '._'])
    return sorted([f for f in os.listdir(path) if f[:1] != '.'])

# take a filename or dirname and return the embedded cutoff value as a float
def dname2cutoff(dname):
    # get stuff before the dash
    cutoff = dname.split('-')[0].replace('_', '.')
    # remove alpha chars
    cutoff = ''.join([char for char in list(cutoff) if (char.isdigit() or char == '.')])
    #print >> sys.stderr, [char for char in list(cutoff) if (char.isdigit() or char == '.')]
    return float(cutoff)

# take name of a PCOC results file, return list of PPs
def pcocResults2PPs(filename):
    results = pd.read_table(filename)
    return list(results["PCOC"])

# take a list and a threshold, return fraction of the list above that threshold
def rateAbove(dataset, threshold):
    binList = [0b1 if i >= threshold else 0b0 for i in dataset]
    return float(np.mean(binList))

# take a list and a threshold, return fraction of the list above that threshold
def rateBelow(dataset, threshold):
    binList = [0b1 if i < threshold else 0b0 for i in dataset]
    return float(np.mean(binList))

# go down the passed PCOC sim-results directory,
# return a DataFrame of positive or negative rates
# (pos=True gives rate above threshold, False gives rate below)
def buildAccuracyTable(dir, threshold=0.8, pos=True):
    accuracyDf = pd.DataFrame() # init DataFrame
    for simDir in listdirSafe(dir):
        workPath = dir + '/' + simDir
        for detDir in listdirSafe(workPath):
            workPath = workPath + '/' + detDir
            # get positive rate of latest run
            latestRun = workPath + '/' + listdirSafe(workPath)[-1]
            resultsFilename = latestRun + '/' + [fname for fname in listdirSafe(latestRun) if fname.find(".results.tsv") > 0][0]
            if pos:
                rate = rateAbove(pcocResults2PPs(resultsFilename), threshold)
            else:
                rate = rateBelow(pcocResults2PPs(resultsFilename), threshold)
            accuracyDf.at[dname2cutoff(simDir), str(dname2cutoff(detDir))] = rate
            # go back up a level
            workPath = dir + '/' + simDir

    return accuracyDf

def main(argv, wayout):
    if not len(argv):
        argv.append('-h')
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    parser.add_argument("-a", "--anc_dir", type=str, help="directory of negative (ancestral) sims", required=True)
    parser.add_argument("-c", "--con_dir", type=str, help="directory of positive (convergent) sims", required=True)
    parser.add_argument("-t", "--threshold", type=float, help="calling threshold for +/-", required=True)
    parser.add_argument("-hm", "--heatmaps", type=str, default=None, help="where to save heatmap PDFs")

    global args
    args = parser.parse_args(argv)

    tpr = buildAccuracyTable(args.con_dir, threshold=args.threshold, pos=True)
    fpr = buildAccuracyTable(args.anc_dir, threshold=args.threshold, pos=True)
    #tnr = buildAccuracyTable(args.anc_dir, threshold=args.threshold, pos=False)
    #fnr = buildAccuracyTable(args.con_dir, threshold=args.threshold, pos=False)

    if args.heatmaps:
        for i, accTable in enumerate((tpr, fpr)):
            # Make a heatmap figure. Graphics parameters are all set in `pcoc_cont_heatmap.py`
            rainbow = True  # Is the trait in question a visible wavelength to be plotted chromatically?
            heatmapper.heatMapDF(accTable, args.heatmaps + "/accuracyHeatmap" + str(i) + ".pdf", rainbow=rainbow)

    return tpr, fpr

if __name__ == "__main__":
    # when called from command line, print results to stdout
    for table in main(sys.argv[1:], sys.stderr):
        print >> sys.stdout, table.to_csv(sep='\t')