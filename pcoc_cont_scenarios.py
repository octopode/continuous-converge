#!/usr/bin/python
#  pcoc_cont_scenarios.py
#
#  This is a wrapper for pcoc_det.py that attempts to adapt the method for continuous traits.
#  Because PCOC requires and old version of Bpp, it is recommended that this be run in an appropriate Docker container.
#
#  If you need to recreate this container:
#  Download the PCOC docker image: docker pull carinerey/pcoc
#  Launch it into bash: docker run -it carinerey/pcoc /bin/bash
#  run pip install pandas biopython matplotlib
#  commit the image, e.g. docker commit 4ff2d5a930dc pcoc-matplot
#
#  example usage:
#  xvfb-run python /data/bin/pcoc/src/pcoc_cont_scenarios.py -t tree/FPs_knowngenes_MAFFT-gblocks-RAxML.tree -o $PWD -c trait/FPemWL_known_noref.tab -aa seq/FPs_knowngenes_MAFFT_339.fasta -d
#
#  Copyright 2018 Jacob Winnikoff <jwinnikoff@mbari.org>
#  uses material by (c) 2017 Carine Rey <carine.rey@ens-lyon.fr>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#

# pip modules
import os
import sys
import argparse
import configparser
import logging
import subprocess
import datetime
import numpy as np
import pandas as pd
from ete3 import Tree, NodeStyle, TreeStyle, TextFace, CircleFace
from pprint import pprint
from Bio import AlignIO

# custom modules
import pcoc_treeviz as treeviz
import pcoc_plots as pviz
from wl2rgb import wavelength_to_rgb

# create logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

### Parse command line args
def parse_args(argv):

    # REQUIRED ARGS #

    parser = argparse.ArgumentParser(prog="pcoc_cont_scenarios.py", description='')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    requiredOptions = parser.add_argument_group('Required arguments')
    requiredOptions.add_argument('-t', "--tree", type=str, help='input tree name', required=True)
    requiredOptions.add_argument('-c', '--cont_trait_table', type=str, help="trait table name", required=True)
    requiredOptions.add_argument('-o', '--output', type=str, help="Output directory", required=True)

    # OPTIONAL ARGS FOR CONTINUOUS PCOC #

    Options = parser.add_argument_group('Options')
    Options.add_argument('-aa', '--aa_align', type=str, help="AA alignment name")
    Options.add_argument('-i', '--invert_trait', action="store_true",
                         help="Invert the binary trait, i.e. assert that low trait value is the convergent state")
    Options.add_argument('-nb', '--num_bins', type=int, default=0,
                         help="Average continuous trait into n equal-width bins, 0 = no binning")
    Options.add_argument('-bw', '--bin_width', type=float, default=0,
                         help="Average continuous trait into bins of width w, 0 = no binning")
    Options.add_argument('--decimal', type=int, default=0,
                         help="Round traits to n decimals in filenames and figures")
    Options.add_argument('-tt', '--test_trees', action="store_true",
                         help="Draw test trees to evaluate the discretization scheme")
    Options.add_argument('--det', action="store_true", help="Set to actually run pcoc_det.py")
    Options.add_argument('--sim', type=float, default=0, help="Set PCOC PP threshold for post hoc simulation, 0.0001 = all sites")
    Options.add_argument('-k', '--key_seq', type=str, help="Name of key sequence on which to index the output columns")
    Options.add_argument('-m', '--master_table', type=str, help="Save collated master data table at...")
    Options.add_argument('-hm', '--heatmap', type=str,
                         help="Render heatmap from the latest available set of data and save it here")
    Options.add_argument('-mp', '--manhattan', type=str, help="Save Manhattan plot at...")
    Options.add_argument('-d', '--debug', action="store_true", help="debug mode")

    # parse the above from command line
    contArgs, passthruArgv = parser.parse_known_args(argv[1:])

    # ARGS FOR SITE DETECTION #

    detParser = argparse.ArgumentParser(prog="pcoc_det.py", description='optional args to pass to pcoc_det.py')
    detOptions = detParser.add_argument_group('Arguments for site detection')

    detOptions.add_argument('-f', '--filter_t', type=float,
                              help="ALL model: Posterior probability threshold to put result in \"filtered\" results. (default: 0.99)",
                              default=-1)
    detOptions.add_argument('-f_pcoc', '--filter_t_pcoc', type=float,
                              help="PCOC model: Posterior probability threshold to put result in \"filtered\" results. If = -1, take the value of -f, if > 1, discard this model. (default: -1)",
                              default=-1)
    detOptions.add_argument('-f_pc', '--filter_t_pc', type=float,
                              help="PC model: Posterior probability threshold to put result in \"filtered\" results. If = -1, take the value of -f, if > 1, discard this model.(default: -1)",
                              default=-1)
    detOptions.add_argument('-f_oc', '--filter_t_oc', type=float,
                              help="OC model: Posterior probability threshold to put result in \"filtered\" results. If = -1, take the value of -f, if > 1, discard this model.(default: -1)",
                              default=-1)
    detOptions.add_argument('-ph', type=str,
                              help="Add these positions in the filtered position and highlight them with a star in the plot",
                              default=False)
    detOptions.add_argument('--plot', action="store_true",
                              help="Plot the tree and the filtered sites of the alignment with their corresponding score.",
                              default=False)
    detOptions.add_argument('--plot_complete_ali', action="store_true",
                              help="Plot the tree and each site of alignment with its corresponding score. (Can take time to be openned)",
                              default=False)
    detOptions.add_argument('-plot_title', type=str,
                              help="Title of each plot (default: None)",
                              default="")
    detOptions.add_argument('--reorder', action="store_true",
                              help="reorder the filtered plot by score.categories (>= 0.99, >=0.9, >= 0.8, < 0.8)",
                              default=False)
    detOptions.add_argument('--svg', action="store_true",
                              help="additional svg output plots.",
                              default=False)
    detOptions.add_argument('--no_cleanup_fasta', action="store_true",
                              help="Do not cleanup the fasta directory after the run.",
                              default=False)
    detOptions.add_argument('-CATX_est', type=int, choices=[10, 60],
                                 help="Profile categorie to estimate data (10->C10 or 60->C60). (default: 10)",
                                 default=60)
    detOptions.add_argument('--gamma', action="store_true",
                                 help="Use rate_distribution=Gamma(n=4) instead of Constant()",
                                 default=False)
    detOptions.add_argument('--inv_gamma', action="store_true",
                                 help="Use rate_distribution=Gamma(n=4) instead of Constant()",
                                 default=False)
    detOptions.add_argument('--max_gap_allowed', type=int,
                                 help="max gap allowed to take into account a site (in %%), must be between 0 and 100 (default:5%%)",
                                 default=5)
    detOptions.add_argument('--max_gap_allowed_in_conv_leaves', type=int,
                                 help="max gap allowed in convergent leaves to take into account a site (in %%), must be between 0 and 100 (default:5%%)",
                                 default=5)
    detOptions.add_argument('--no_cleanup', action="store_true",
                                 help="Do not cleanup the working directory after the run.",
                                 default=False)
    detOptions.add_argument("-LD_LIB", metavar='LD_LIBRARY_PATH', type=str, default="",
                                 help="Redefine the LD_LIBRARY_PATH env variable, bppsuite library must be present in the $PATH and in the LD_LIBRARY_PATH")

    # parse the above and leave remaining argv for sim
    detArgs, simArgv = detParser.parse_known_args(passthruArgv)

    # ARGS FOR SIMULATION #

    simParser = argparse.ArgumentParser(prog="pcoc_sim.py", description='optional args to pass to pcoc_sim.py')
    simOptions = simParser.add_argument_group('Arguments for simulation-based error control')

    simOptions.add_argument('-nb_sampled_couple', type=int, metavar="INT",
                             help="For each convergent scenario, number of simulated alignment with different sampled couple of profiles (Ancestral/Convergent). (default: 1)",
                             default=1)
    simOptions.add_argument('-n_sites', type=int, metavar="INT",
                             help="Number of simulated sites per alignment. (default: 1000)",
                             default=1000)
    simOptions.add_argument('-CATX_sim', type=int, choices=[10, 60],
                             help="Profile categories to simulate data (10->C10 or 60->C60). (default: 60)",
                             default=60)
    simOptions.add_argument('-min_dist_CAT', type=float, metavar="FLOAT",
                             help="Minimum distance between Ancestral and Convergent profiles to simulate the alignment (default: no limits)",
                             default=0)
    simOptions.add_argument('--plot_ali', action="store_true",
                             help="For each couple of profiles, plot a summary of the convergent scenario containing the tree and the alignment.",
                             default=False)
    simOptions.add_argument('--get_likelihood_summaries', action="store_true",
                             help="For each couple of profiles, write a summary of the likelihoods per site.",
                             default=False)
    simOptions.add_argument('--no_clean_seqs', action="store_true",
                             help="Do not cleanup the sequences after the run.",
                             default=False)

    # parse the above and leave remaining argv for det
    simArgs, detArgv = simParser.parse_known_args(passthruArgv)

    return contArgs, detArgv, simArgv


def main(contArgs, detArgv, simArgv):
    ##########
    # inputs #
    ##########
    startDateTime = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

    ### Set up the log file
    LogFile = contArgs.output + "/det_{}.log".format(startDateTime)

    ### Set up the loggers
    det_logger = logging.getLogger("pcoc_det_mod")

    # create file handler which logs even debug messages
    fh = logging.FileHandler(LogFile)
    fh.setLevel(logging.DEBUG)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    if contArgs.debug:
        ch.setLevel(logging.DEBUG)
    else:
        ch.setLevel(logging.INFO)
    # create formatter and add it to the handlers
    formatter_fh = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    formatter_ch = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter_fh)
    ch.setFormatter(formatter_ch)
    # add the handlers to the loggers
    logger.addHandler(fh)
    logger.addHandler(ch)
    det_logger.addHandler(fh)
    det_logger.addHandler(ch)
    # log the calling command
    logger.debug(sys.argv)

    # Check treefile
    if not os.path.isfile(contArgs.tree):
        #print >> sys.stderr, "# ERROR: {} does not exist".format(contArgs.tree)
        logger.error("{} does not exist".format(contArgs.tree))
        sys.exit(1)

    # Load treefile
    tree = init_tree(contArgs.tree)
    # not mucking with additive trees yet; ultrametrize the tree and normalize to length 1
    tree.convert_to_ultrametric(tree_length=1)

    ### Generate convergent scenarios

    # Load in dict of traits keyed on species. Note hardcoding of 'sp' colName!
    tipTraits = pd.read_table(contArgs.cont_trait_table).set_index('sp').transpose().to_dict(orient='r')[0]
    # Do BM trait reconstruction
    nodeTraits = ancR(tree, tipTraits)
    # get cutoffs using specified binning
    cutoffs = binBy(nodeTraits, fixNumBins=contArgs.num_bins, fixBinWidth=contArgs.bin_width)[0]
    # Generate discrete trait matrix for all cutoff values
    binDF = discretize(nodeTraits, cutoffs, contArgs.invert_trait)
    # consolidate redundant scenarios
    # binDF = uniqueColumnFilter(binDF)[0] # element [1] is the chaff
    binDF = uniqueScenarios(binDF)
    # eliminate convergent root scenarios
    # binDF = convergentRootFilter(binDF)[0] # element [1] is the chaff
    binDF = convergentRootFilter(binDF)
    # convert binary dataframe into scenario strings
    scenarios = scenarioStrings(binDF, tree)

    scenDir = contArgs.output + "/Scenarios"

    ### Draw trees depicting convergent scenarios

    if contArgs.test_trees or contArgs.det:

        # make a dir for the test trees
        testDir = contArgs.output + "/TestTrees"
        try:
            os.mkdir(testDir)
        except:
            pass

        #print >> sys.stderr, "# MESSAGE: Drawing trees..."
        logger.info("Drawing trees...")
        # draw trait-colored tree
        treeviz.traitTree(tree, nodeTraits, wavelength_to_rgb, testDir)
        # treeviz.traitTreeMinimal(tree, nodeTraits, wavelength_to_rgb, contArgs.output)  # draw trait-colored tree
        # draw boring test trees that indicate convergent scenario with red/orange
        treeviz.testTrees(tree, scenarios, testDir, contArgs.tree, contArgs.decimal)
        # treeviz.testTreesMinimal(tree, scenarios, nodeTraits, wavelength_to_rgb, contArgs.tree, contArgs.output) # draw test trees with color-coded nodes

    ### run site detection

    if contArgs.det:
        # make a dir for all that pcoc_det output
        try:
            os.mkdir(scenDir)
        except:
            pass

        # log the continuous trait values and convergent scenarios used
        logger.debug("continuous trait values:\n{}".format(tipTraits))

        # run pcoc_det.py for each convergent scenario
        detArgvStatic = ["pcoc_det.py",
                   "-t", contArgs.tree,
                   "-aa", contArgs.aa_align,
                   "--no_cleanup"] + detArgv # add additional site-detect args from command line
        # error will be thrown if any of these are specified redundantly at the command line

        # DETECT LOOP #
        logger.info("Detecting convergent sites...")
        for num, cutoff in enumerate(sorted(scenarios.keys())):
            scenarioString = scenarios[cutoff]
            logger.info("convergent scenario {}/{}: cutoff {}: {}".format(num + 1, len(scenarios), cutoff, scenarioString))
            detArgvDynamic = detArgvStatic + [
                "-o", scenDir + '/' + str(cutoff).replace('.', '_'),
                "-m", str(scenarios[cutoff])]

            # run pcoc.det.py
            subprocess.call(detArgvDynamic)
        # END DETECT LOOP #

    ### Collate, visualize results

    # assemble long-format master DataFrame
    # columns:
    # Site
    # PCOC
    # PC
    # OC
    # Cutoff
    # Scenario
    # Profile_A
    # Profile_C
    # NB_CAT

    metaDf = consolidatePCOCOutput(scenarios, scenDir, contArgs.aa_align)

    # map the PP scores to key sequence positions if desired
    # alignment reindexing is done post-analysis so user can switch key seqs on the fly without reanalyzing
    #TODO: revise key seq mapping so that all renumbered files are saved as copies and original is always preserved
    '''
    if contArgs.key_seq:
        alignment = AlignIO.read(contArgs.aa_align, format="fasta")
        # get key columns of the alignment
        keyCols = keyAlignmentColumns(alignment, contArgs.key_seq)
        # filter the 2-column data table to contain only these sites
        metaDf = metaDf[metaDf["Sites"].isin(keyCols)]
        # reindex the sites
        metaDf["Sites"] = [keyCols.index(i) + 1 for i in metaDf["Sites"]]  # note index shift!
    '''

    # recast the molten/long-format metaDF to wide format
    metaDfWide = metaDf.pivot(columns="Cutoff", index="Sites", values="PCOC")
    # the transpose of this is used for plotting
    heatmapDf = metaDf.pivot(index="Cutoff", columns="Sites", values="PCOC")

    # contArgs.sim is a float value specifying how strong the a priori PP signal needs to be at a site to warrant
    # post hoc simulation.
    if contArgs.sim:
        # init post hoc simulation DF with sites index
        simDf = pd.DataFrame(index = metaDfWide.index)

        # get the highest PP for each site and the cutoff for that PP
        simDf["PP_max"] = [max(row[1].tolist()) for row in metaDfWide.iterrows()]

        # get the cutoff for that PP
        simDf["Cutoff"] = [metaDfWide.columns[row[1].tolist().index(simDf["PP_max"][rowNum+1])] for rowNum, row in enumerate(metaDfWide.iterrows())]
        # below should be the same command, using .iloc() for speed
        #simDf["Cutoff"] = [metaDfWide.columns[row[1].tolist().index(simDf.iloc(rowNum + 1, "PP_max"))] for rowNum, row in enumerate(metaDfWide.iterrows())]

        # lnL matrix lookup is not by site, but by cutoff
        # so lnL matrices for all top cutoffs are loaded into memory
        # if this gets too heavy I can save/read the matrices to disk
        logger.info("Loading site-model log-likelihood matrices...")
        # make a set of just the relevant cutoffs
        maxPPcutoffsSet = set([cutoff for cutoff in simDf["Cutoff"] if cutoff >= contArgs.sim])
        # load 3D matrices for those cutoffs
        # this is a slow step, probably from all the file handling
        # its complexity should scale with number of cutoffs, not sites.
        lnLMatrices = {cutoff: consolidateBPPOutput(cutoff, scenDir) for cutoff in maxPPcutoffsSet}

        # iterate over the sites and assign profiles and model lnLs to each
        logger.info("Recovering ML CATegories for sites...")
        for site in range(simDf.shape[0]):
            # if the site made the cut for simulation
            if simDf["Cutoff"].tolist()[site] >= contArgs.sim:
                # get the profiles and the model lnL and tack them on the row
                simDf.loc[site + 1, "CAT_Anc"], simDf.loc[str(site + 1), "CAT_Con"], simDf.loc[site + 1, "Model_lnL"] = getMLCATProfiles(site, lnLMatrices[cutoff])

        postHocFilename = os.path.splitext(contArgs.master_table)[0] + "_posthoc.tsv"
        simDf.to_csv(postHocFilename, sep='\t')

    if contArgs.heatmap:
        #print >> sys.stderr, "# MESSAGE: Drawing heatmap:"
        # Make a heatmap figure. Graphics parameters are all set in `pcoc_cont_heatmap.py`
        heatmapPath = contArgs.output + '/' + contArgs.heatmap
        rainbow = True  # Is the trait in question a visible wavelength to be plotted chromatically?
        pviz.heatMapDF(heatmapDf, heatmapPath, rainbow=rainbow)
        # tell user where it was put
        #print >> sys.stderr, heatmapPath
        logger.info("Drawing heatmap at {}".format(heatmapPath))

    if contArgs.manhattan:

        ### Sum columns of non-exclusive probabilities
        def sumProbs(probArray):
            # the array casting is important, because numpy uses the operators element-wise
            sums = np.array([0.] * probArray.shape[1])
            for row in probArray:
                sums = sums + np.array(row) - (sums * np.array(row))
            return sums

        #print >> sys.stderr, "# MESSAGE: Drawing Manhattan plot:"
        manhattanSeries = sumProbs(np.array(heatmapDf))  # convert to array
        manhattanPath = contArgs.output + '/' + contArgs.manhattan
        pviz.manhattanPlot(manhattanSeries, manhattanPath, contArgs.key_seq)
        # tell user where it was put
        #print >> sys.stderr, manhattanPath
        logger.info("Drawing Manhattan plot at {}".format(manhattanPath))

    if contArgs.master_table:
        #print >> sys.stderr, "# MESSAGE: Saving master table of results:"
        # Print master data table. Note that it is transpose of the table sent to heatMapDF()
        metaDfWide.to_csv(contArgs.master_table, sep='\t')
        # tell user where it was put
        #print >> sys.stderr, contArgs.master_table
        logger.info("Saving master table of results at {}".format(contArgs.master_table))

## HELPER FUNCTIONS ##

# take an argparse.Namespace and convert it back to an argv-style list
def args2argv(args):

    # get (flag, arg) tuples
    argTupls = args.__dict__.items()
    # prepend '-- ' to the flags
    argTupls = [('--'+el[0], el[1]) for el in argTupls]
    # return flattened argv list
    return [str(el) for tupl in argTupls for el in tupl]

### Perform BM continuous-trait ancestral state reconstruction. Return list of trait values indexed on node #.
def ancR(tree, tipTraits):

    # this command does something funky to the branch lengths!
    #tree.convert_to_ultrametric(tree_length=1)

    # Check that the trait list has values for all leaves
    dataError = False
    for leaf in tree.get_leaves():
        if leaf.name not in tipTraits.keys():
            #print >> sys.stderr, "No trait data for {}!".format(leaf.name)
            logger.error("No trait data for {}!".format(leaf.name))
            # makes sure all the missing taxa get read thru before the program crashes
            dataError = True
    if dataError:
        #print >> sys.stderr, "exiting"
        sys.exit(1)

    # Init list with an element for each node in tree
    nodeTraits = [None] * len(list(tree.traverse()))
    nodeId = 0
    for node in tree.traverse("postorder"): # Amend this to use enumerate()

        if node.name in tipTraits.keys(): # if the node is a tip or its state is specified
            fixTrait = tipTraits[node.name]
            node.add_features(trait = fixTrait)
            nodeTraits[nodeId] = fixTrait

        else: # if the state needs to be reconstructed
            daughters = node.get_descendants()[:2]
            daughterTimes = [tree.get_distance(node, daughter) for daughter in daughters]
            daughterTraits = [daughter.trait for daughter in daughters]
            rate = (daughterTraits[1] - daughterTraits[0]) / sum(daughterTimes)
            recTrait = daughterTraits[0] + (rate * daughterTimes[0])

            node.add_features(trait = recTrait)
            nodeTraits[nodeId] = recTrait

        nodeId += 1

    return nodeTraits


# take a list of node traits, bin them as specified, return both inter-bin cutoffs and bin means
def binBy(nodeTraits, fixNumBins = 0, fixBinWidth = 0):
    #TODO: implement unbiased rounding for np.digitize()

    # can only specify number of bins OR bin width
    if fixNumBins and fixBinWidth:
        sys.exit(1)

    # with specified number of bins
    if fixNumBins:
        numBins = fixNumBins
        binWidth = float(max(nodeTraits) - min(nodeTraits)) / fixNumBins
        binBounds = [min(nodeTraits) + (binWidth * i) for i in range(fixNumBins + 1)]
        binBounds[-1] *= 1.001 # include the right boundary
        bIndices = [i-1 for i in np.digitize(nodeTraits, binBounds)]
        activeBins = sorted(list(set(bIndices)))
        #print >> sys.stderr, "# MESSAGE: Use of {} equal-width bins yielded {} unique trait values:".format(fixNumBins, len(activeBins))
        logger.info("Use of {} equal-width bins yielded {} unique trait values:".format(fixNumBins, len(activeBins)))

    # with specified bin width, centered on midpoint of range
    if fixBinWidth:
        numBins = int((max(nodeTraits) - min(nodeTraits)) / fixBinWidth) + 1
        midpoint = float(max(nodeTraits) + min(nodeTraits)) / 2.0
        start = midpoint - (float(fixBinWidth) * numBins / 2.0)
        binBounds = [start + (fixBinWidth * i) for i in range(numBins + 1)]
        binBounds[-1] *= 1.001  # include the right boundary
        bIndices = [i-1 for i in np.digitize(nodeTraits, binBounds)]
        activeBins = sorted(list(set(bIndices)))
        #print >> sys.stderr, "# MESSAGE: Use of bins {} trait units wide yielded {} unique trait values:".format(fixBinWidth, len(activeBins))
        logger.info("Use of bins {} trait units wide yielded {} unique trait values:".format(fixBinWidth, len(activeBins)))

    # now, bin the values and average them
    if fixNumBins or fixBinWidth:
        toAverage = [list()] * numBins
        for bin, val in zip(bIndices, nodeTraits):
            toAverage[bin] = toAverage[bin] + [val]
        # only keep the bins with stuff in them
        toAverage = [bin for bin in toAverage if bin]

        # place cutoffs in between the bin boundaries
        cutoffs = [np.mean((max(toAverage[i]), min(toAverage[i + 1]))) for i in range(len(toAverage) - 1)]

        # write bin means
        nodeTraits = [np.mean(bin) for bin in toAverage]

        # report on binning operation
        for avg, init in zip(nodeTraits, toAverage):
            #print >> sys.stderr, "# {} <- {}".format(avg, init)
            logger.info("{} <- {}".format(avg, init))
    # if no binning was specified, place cutoffs between unique trait values
    else:
        nodeTraits = np.unique(nodeTraits)
        cutoffs = [np.mean((nodeTraits[i], nodeTraits[i+1])) for i in range(len(nodeTraits)-1)]

    return cutoffs, nodeTraits


# Take a list of node traits and return a boolean DataFrame with all possible discretizations of the list.
# DF headers are the discretization cutoffs.
def discretize(nodeTraits, cutoffs, invert = False):

    # get unique node trait values and sort ascending
    uniqueNodeTraits = sorted(np.unique(nodeTraits))

    # get all possible cutoffs
    #cutoffs = [np.mean((uniqueNodeTraits[i], uniqueNodeTraits[i+1])) for i in range(len(uniqueNodeTraits)-1)]

    # init pd.DataFrame
    binDF = pd.DataFrame()
    # put continuous values into the first column. Rows are nodes in post-order.
    binDF["trait_cont"] = nodeTraits

    # append binary columns to DataFrame
    for cutoff in cutoffs:
        binDF[cutoff] = pd.cut(binDF["trait_cont"], [min(nodeTraits)]+[cutoff]+[max(nodeTraits)],
                                    include_lowest=True, labels=False)
    # then drop the first column
    binDF = binDF.drop(binDF.columns[0], axis=1)

    # convert binary values to Boolean
    if invert:
        # invert the trait if specified
        binDF = ~binDF.astype(bool)
    else:
        # otherwise high trait value is considered the convergent state
        binDF = binDF.astype(bool)

    return binDF


# Take a discrete trait DataFrame, consolidate indentical columns and average the headers (which should be floats)
def uniqueScenarios(binDF):

    # make a set of (unique) column tuples
    uniqueCols = list(set([tuple(binDF[column]) for column in binDF.columns]))

    numInitCols = len(binDF.columns)
    numUniqueCols = len(uniqueCols)

    if numUniqueCols < numInitCols:

        # group the headers matching each unique column in a list of lists
        toAverage = [list()] * len(uniqueCols)
        for i, col in enumerate(uniqueCols):
            for colName, series in sorted(list(binDF.iteritems())):
                if tuple(series) == col:
                    # append-in-place no good here
                    toAverage[i] = toAverage[i] + [colName]

        # average each list in the list of lists
        avgCutoffs = [np.mean(cuts) for cuts in toAverage]

        # list them
        #print >> sys.stderr, "# MESSAGE: {} scenarios were consolidated into {} unique scenarios:".format(numInitCols,
        #                                                                                                  numUniqueCols)
        logger.info("{} scenarios were consolidated into {} unique scenarios:".format(numInitCols, numUniqueCols))
        for avg, init in zip(avgCutoffs, toAverage):
            #print >> sys.stderr, "# {} <- {}".format(avg, init)
            logger.info("{} <- {}".format(avg, init))

        # construct consolidated dataframe
        binDF = pd.DataFrame(columns = avgCutoffs, data = zip(*uniqueCols))
        # sort the new averaged columns
        binDF.sort_index(axis = 1, inplace = True)

    else:

        #print >> sys.stderr, "# MESSAGE: Scenarios for all cutoffs are unique"
        logger.info("Scenarios for all cutoffs are unique.")

    return binDF


### Remove and report on scenarios in which the root is called convergent
def convergentRootFilter(binDF):

    # store the original columns
    origCutoffs = binDF.columns
    # get the columns whose last (root) element is True
    convergentCutoffs = [colName for colName in origCutoffs if binDF[colName].iloc[-1]]
    # drop those columns
    binDF.drop(labels = convergentCutoffs, axis = 1, inplace = True)

    # notify user of outcome of root filter
    #print >> sys.stderr, "# MESSAGE: {}/{} scenarios eliminated due to convergent root:".format(len(convergentCutoffs), len(origCutoffs))
    logger.info("{}/{} scenarios eliminated due to convergent root:\n{}".format(len(convergentCutoffs), len(origCutoffs), convergentCutoffs))

    if len(convergentCutoffs) >= len(origCutoffs) / 2:
        #print >> sys.stderr, "# WARNING: Consider inverting your trait!"
        logger.warning("Consider inverting your trait!")

    return binDF


### Use tree and binary trait matrix to generate 'convergent scenario' strings readable by PCOC scripts
def scenarioStrings(binDF, tree):

    # scenario strings keyed on cutoff value
    scenarios = dict()

    for col in binDF.columns:
        scenario = str()
        for node, state in enumerate(binDF[col]):
            # if convergent state is true
            if state:
                # if not the root node
                if (node < binDF.shape[0] - 1):
                    parentND = tree.search_nodes(ND=node)[0].up.ND
                    if binDF[col][parentND]:
                        scenario = ',' + str(node) + scenario
                    else:
                        scenario = '/' + str(node) + scenario
                else:
                    scenario = '/' + str(node) + scenario

        # remove leading '/' and add scenario to dict, keyed on cutoff
        scenarios[col] = scenario[1:]

    # return dict with cutoff: scenarioString
    return scenarios

def scenarioString2Dict(scenarioString):
    manual_mode_nodes = {"T": [], "C": []}
    p_events = scenarioString.strip().split("/")
    for e in p_events:
        l_e = map(int, e.split(","))
        manual_mode_nodes["T"].append(l_e[0])
        manual_mode_nodes["C"].extend(l_e[1:])

    return p_events, manual_mode_nodes

### Gather up all the most recent output files in ScenDir and put them into a master dataframe
def consolidatePCOCOutput(scenarios, scenDir, aa_align):

    # init output DF
    metaDf = pd.DataFrame()

    # get the output file from the latest run
    ppFileBasename = os.path.splitext(os.path.basename(aa_align))[:-1][0] + ".results.tsv"

    for cutoff in sorted(scenarios.keys()):

        try:
            subDir = scenDir + '/' + str(cutoff).replace('.','_')
            latestRun = sorted(os.listdir(subDir))[-1]
            ppFile = subDir + "/" + latestRun + "/" + ppFileBasename
            cutoffPPs = pd.read_table(ppFile)
            # store the cutoff and the scenario string
            cutoffPPs["Cutoff"] = cutoff
            cutoffPPs["Scenario"] = scenarios[cutoff]
            metaDf = metaDf.append(cutoffPPs)
        except:
            logger.warning("Posterior probabilities not loaded for cutoff {}.".format(cutoff))

    return metaDf

### Gather up the .infos files for a particular cutoff in ScenDir and put the lnL data into a 3D array
def consolidateBPPOutput(cutoff, scenDir, withOneChange = True):

    # search criterion for files to load from the Estimations dir
    if withOneChange:
        searchPattern = "withOneChange.infos"
    else:
        searchPattern = "noOneChange.infos"

    try:
        subDir = scenDir + '/' + str(cutoff).replace('.', '_')
        # run folders are named by numeric datetime, so last is most recent
        latestRun = sorted(os.listdir(subDir))[-1]
        estimsDir = subDir + "/" + latestRun + "/Estimations"
        # get list of target filenames
        lnLfnames = [fname for fname in os.listdir(estimsDir) if fname.endswith(searchPattern)]
    except:
        logger.warning("Log likelihoods not loaded for cutoff {}.".format(cutoff))
        return None

    ## build site x aProfile x cProfile array

    # get the profiles and profile counts
    profiles = [infosFname2Profiles(fname) for fname in lnLfnames]
    profilesA = zip(*profiles)[0]
    numProfilesA = max(profilesA)
    profilesC = zip(*profiles)[1]
    numProfilesC = max(profilesC)

    # taste the first file to get number of sites
    numSites = pd.read_table(estimsDir + "/" + lnLfnames[0]).shape[0]

    # init the array
    lnLs = np.empty((numSites, numProfilesA, numProfilesC))

    # iteratively build the array
    for i, fname in enumerate(lnLfnames):
        #print pd.read_table(estimsDir + "/" + lnLfnames[i])["lnL"].tolist()
        lnLs[:][profilesA[i]-1][profilesC[i]-1] = pd.read_table(estimsDir + "/" + lnLfnames[i])["lnL"].tolist()

    return lnLs

# get the ancestral and convergent CAT profiles from name of a PCOC/BPP .infos file
def infosFname2Profiles(fname):
    # split filename on '_'; take the 2nd- and 3rd-to-last elements
    return [int(cat) for cat in fname.split('_')[-3:-1]]

### Take a site (0-indexed) and a np.array of log-likelihoods for different ancestral/convergent CAT models,
### return the ancestral and convergent CATegory pair with the highest log likelihood. Also return lnL for diagnostics.
def getMLCATProfiles(site, lnLs):

    # get A, C indices of the max lnL for a given site in the matrix
    profileA, profileC = np.unravel_index(np.argmax(lnLs[site][:][:]), (lnLs.shape[1], lnLs.shape[2]))
    # get the actual lnL
    lnL = lnLs[site][profileA][profileC]

    # return 1-indexed CAT numbers
    return int(profileA + 1), int(profileC + 1), lnL

### return list of columns of an alignment that are not gaps in the key sequence
def keyAlignmentColumns(algt, keySeqID):

    algtLength = algt.get_alignment_length()

    # Get the key sequence
    for record in algt:
        if record.id == keySeqID:
            keySeq = record.seq

    mappedCols = [None] * len(str(keySeq).replace('-', ''))
    keyIndex = 0
    for i in range(algtLength):  # only iterate over columns that are not gaps in target seq
        if keySeq[i] != "-":  # meaning anything except gaps
            mappedCols[keyIndex] = i + 1  # map the alignment column to the key column. Gotta shift the index!
            keyIndex += 1

    return mappedCols


### Number tree nodes for consistent reference
def init_tree(nf):
    t = Tree(nf)

    for i, n in enumerate(t.traverse("postorder")):
        n.add_features(ND = i)

    return t

if __name__ == "__main__":
    contArgs, detArgv, simArgv = parse_args(sys.argv)
    main(contArgs, detArgv, simArgv)