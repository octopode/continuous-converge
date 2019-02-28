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
# changelog
#
# v1.9 - 20190212: Reimplemented cutoff placement using uniform trait bins to mitigate sampling bias
# v2.0 - 20190212: Eliminated subprocess calls from the script; PCOC routines are now imported as modules
#

import os
import sys
import logging
import subprocess
import pcoc_det_mod as det
#import pcoc_sim_mod as sim
import argparse
import configparser
import datetime
import numpy as np
import scipy.stats as stat
import random
import pandas as pd
from ete3 import Tree, NodeStyle, TreeStyle, TextFace, CircleFace
from pprint import pformat
from Bio import AlignIO
from ast import literal_eval
from wl2rgb import wavelength_to_rgb
import pcoc_treeviz as viz
import pcoc_cont_heatmap as heatmapper

# create logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

def parse_args():
    # parse command line arguments

    ### Option defining
    parser = argparse.ArgumentParser(prog="pcoc_cont_scenarios.py", description='')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    ##############
    requiredOptions = parser.add_argument_group('Required arguments')
    requiredOptions.add_argument('-t', "--tree", type=str, help='input tree name', required=True)
    requiredOptions.add_argument('-o', '--output', type=str, help="Output directory", required=True)
    requiredOptions.add_argument('-c', '--cont_trait_table', type=str, help="trait table name", required=True)
    ##############
    Options = parser.add_argument_group('Options')
    Options.add_argument('-aa', '--aa_align', type=str, help="AA alignment name")
    Options.add_argument('-k', '--key_seq', type=str, help="Name of key sequence on which to index the output columns")
    Options.add_argument('-tt', '--test_trees', action="store_true",
                         help="Draw test trees to evaluate the discretization scheme")
    Options.add_argument('-i', '--invert_trait', action="store_true",
                         help="Invert the binary trait, i.e. assert that low trait value is the convergent state")
    Options.add_argument('-d', '--det', action="store_true", help="Set to run site detection using passed alignment")
    Options.add_argument('-s', '--sim', type=str, default=None, help="Location of error curves for simulated alignments")
    Options.add_argument('-f', '--float', type=int, default=0,
                         help="Round traits to f decimals in filenames and figures")
    Options.add_argument('-nb', '--num_bins', type=int, default=0,
                         help="Average continuous trait into n equal-width bins, 0 = no binning")
    Options.add_argument('-bw', '--bin_width', type=float, default=0,
                         help="Average continuous trait into bins of width w, 0 = no binning")
    # Options.add_argument('-p', '--precision', type=float, default=0.0,
    #                     help="Minimum difference in trait cutoffs between consecutive scenarios to be tested")
    Options.add_argument('-hm', '--heatmap', type=str,
                         help="Render heatmap from the latest available set of data and save it here")
    Options.add_argument('-m', '--master_table', type=str, help="Save collated master data table at...")
    # TODO: have master table go to stdout, but for some reason my stderr goes there too!
    Options.add_argument('-mp', '--manhattan', type=str, help="Save Manhattan plot at...")
    ##############

    ##############
    AdvancedOptions = parser.add_argument_group('OPTIONS FOR ADVANCED USAGE')
    AdvancedOptions.add_argument('-CATX_est', type=int, choices=[10, 60],
                                 help="Profile categorie to estimate data (10->C10 or 60->C60). (default: 10)",
                                 default=10)
    AdvancedOptions.add_argument('--gamma', action="store_true",
                                 help="Use rate_distribution=Gamma(n=4) instead of Constant()",
                                 default=False)
    AdvancedOptions.add_argument('--inv_gamma', action="store_true",
                                 help="Use rate_distribution=Gamma(n=4) instead of Constant()",
                                 default=False)
    AdvancedOptions.add_argument('--max_gap_allowed', type=int,
                                 help="max gap allowed to take into account a site (in %%), must be between 0 and 100 (default:5%%)",
                                 default=5)
    AdvancedOptions.add_argument('--max_gap_allowed_in_conv_leaves', type=int,
                                 help="max gap allowed in convergent leaves to take into account a site (in %%), must be between 0 and 100 (default:5%%)",
                                 default=5)
    AdvancedOptions.add_argument('--no_cleanup', action="store_true",
                                 help="Do not cleanup the working directory after the run.",
                                 default=False)
    AdvancedOptions.add_argument('--no_cleanup_fasta', action="store_true",
                              help="Do not cleanup the fasta directory after the run.",
                              default=False)
    AdvancedOptions.add_argument("-LD_LIB", metavar='LD_LIBRARY_PATH', type=str, default="",
                                 help="Redefine the LD_LIBRARY_PATH env variable, bppsuite library must be present in the $PATH and in the LD_LIBRARY_PATH")
    AdvancedOptions.add_argument('--debug', action="store_true",
                                 help="debug mode",
                                 default=False)

    ##############
    Options_ben = parser.add_argument_group('Detection options')


    ### Option parsing
    return parser.parse_args()

# A separate argparser to handle simulation parameters
def parse_args_sim():

    ### Option defining
    parser = argparse.ArgumentParser(prog="pcoc_cont_scenarios.py", description='')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    ##############
    Options_ali = parser.add_argument_group('Alignment simulation options')
    Options_ali.add_argument('-nb_sampled_couple', type=int, metavar="INT",
                             help="For each convergent scenario, number of simulated alignment with different sampled couple of profiles (Ancestral/Convergent). (default: 1)",
                             default=1)
    Options_ali.add_argument('-n_sites', type=int, metavar="INT",
                             help="Number of simulated sites per alignment. (default: 100)",
                             default=100)
    Options_ali.add_argument('-CATX_sim', type=int, choices=[10, 60],
                             help="Profile categories to simulate data (10->C10 or 60->C60). (default: 60)",
                             default=60)
    Options_ali.add_argument('-min_dist_CAT', type=float, metavar="FLOAT",
                             help="Minimum distance between Ancestral and Convergent profiles to simulate the alignment (default: no limits)",
                             default=0)
    Options_ali.add_argument('--plot_ali', action="store_true",
                             help="For each couple of profiles, plot a summary of the convergent scenario containing the tree and the alignment.",
                             default=False)
    Options_ali.add_argument('--get_likelihood_summaries', action="store_true",
                             help="For each couple of profiles, write a summary of the likelihoods per site.",
                             default=False)
    Options_ali.add_argument('--no_clean_seqs', action="store_true",
                             help="Do not cleanup the sequences after the run.",
                             default=False)
    ##############

    ### Option parsing
    return parser.parse_known_args()


def main():

    ##########
    # inputs #
    ##########
    startDateTime = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    date = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

    args = parse_args()
    if args.sim: args_sim = parse_args_sim()

    ### Set up the log file
    LogFile = args.output + "/det_{}.log".format(startDateTime)

    ### Set up the loggers
    det_logger = logging.getLogger("pcoc_det_mod")
    bpp_logger = logging.getLogger("pcoc.bpp_lib")

    # create file handler which logs even debug messages
    fh = logging.FileHandler(LogFile)
    fh.setLevel(logging.DEBUG)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    if args.debug:
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
    bpp_logger.addHandler(fh)
    bpp_logger.addHandler(ch)
    # log the calling command
    logger.debug(sys.argv)

    # Check treefile
    if not os.path.isfile(args.tree):
        logger.error("{} does not exist".format(args.tree))
        sys.exit(1)

    # Load treefile
    tree = init_tree(args.tree)
    # not mucking with additive trees yet; ultrametrize the tree and normalize to length 1
    #tree.convert_to_ultrametric(tree_length=1)
    # save the Newick string to RAM
    tree_newick = tree.write()

    # Check MSA
    if not os.path.isfile(args.aa_align):
        logger.error("{} does not exist".format(args.aa_align))
        sys.exit(1)
    # Load MSA
    alignment = AlignIO.read(args.aa_align, format="fasta")

    if args.cont_trait_table is not None:

        ## GENERATE BINARY CONVERGENT SCENARIOS ##

        # Load in dict of traits keyed on species. Note hardcoding of 'sp' colName!
        tipTraits = pd.read_table(args.cont_trait_table).set_index('sp').transpose().to_dict(orient='r')[0]
        # Do BM trait reconstruction
        nodeTraits = ancR(tree, tipTraits, logger)
        # get cutoffs using specified binning
        cutoffs = binBy(nodeTraits, logger, fixNumBins=args.num_bins, fixBinWidth=args.bin_width)[0]
        # Generate discrete trait matrix for all cutoff values
        binDF = discretize(nodeTraits, cutoffs, args.invert_trait)
        # consolidate redundant scenarios
        # binDF = uniqueColumnFilter(binDF)[0] # element [1] is the chaff
        binDF = uniqueScenarios(binDF, logger)
        # eliminate convergent root scenarios
        # binDF = convergentRootFilter(binDF)[0] # element [1] is the chaff
        binDF = convergentRootFilter(binDF, logger)
        # convert binary dataframe into scenario strings
        #scenarios = scenarioStrings(binDF, tree)
        scenarios = scenarioDicts(binDF, tree)

        scenDir = args.output + "/scenarios"

        if args.test_trees:

            ## DRAW TEST TREES ##

            # make a dir for the test trees
            testDir = args.output + "/test_trees"
            try:
                os.mkdir(testDir)
            except:
                pass

            # print >> sys.stderr, "# MESSAGE: Drawing trees..."
            logger.info("Drawing trees...")
            # sys.exit(1)
            # draw trait-colored tree
            viz.traitTree(tree, nodeTraits, wavelength_to_rgb, testDir)
            # traitTreeMinimal(tree, nodeTraits, wavelength_to_rgb, args.output)  # draw trait-colored tree
            # draw boring test trees that indicate convergent scenario with red/orange
            viz.testTrees(tree, scenarios, testDir, args.tree, args.float)
            # testTreesMinimal(tree, scenarios, nodeTraits, wavelength_to_rgb, args.tree, args.output) # draw test trees with color-coded nodes

        if args.det:

            ## DETECT CONVERGENT SITES ##

            # log the continuous trait values used
            #detLogHandle.write("continuous trait values:\n")
            logger.debug("continuous trait values:\n" + pformat(tipTraits))
            #pprint(tipTraits, stream=detLogHandle)

            # get temp directories for pcoc_det
            tempDirs = det.get_temp_dirs(args.aa_align, args.output)

            # init dict of master dataframes
            masters = dict()
            models = ["PCOC", "PC", "OC"]
            for model in models:
                masters[model] = pd.DataFrame()

            # run pcoc_det.py for each convergent scenario
            #print >> sys.stderr, "# MESSAGE: Detecting convergent sites for scenarios..."
            logger.info("Detecting convergent sites for scenarios...")
            # init list of g_tree objects
            g_trees = []
            for num, cutoff in enumerate(sorted(scenarios.keys())):

                logger.info("Processing convergent scenario {}/{}".format(num + 1, len(scenarios)) + ", trait cutoff = {}\n".format(cutoff))

                # run the detect routine
                g_tree, det_df = det.mk_detect(scenarios[cutoff], tree_newick, alignment,
                          *tempDirs,
                          NbCat_Est = args.CATX_est,
                          date = date,
                          max_gap_allowed = args.max_gap_allowed,
                          gamma = args.gamma,
                          inv_gamma = args.inv_gamma
                          )

                # store the g_tree
                g_trees.append(g_tree)

                # put results into the master dataframes
                for model in models:
                    masters[model][cutoff] = det_df[model]

            # once all the cutoffs have been run, clean up the temp dirs
            if not args.no_cleanup:
                det.remove_folder(g_tree.repest)
                det.remove_folder(g_tree.repbppconfig)
                det.remove_folder(g_tree.reptree)
                if not args.no_cleanup_fasta:
                    det.remove_folder(g_tree.repfasta)

            #TEST
            print masters

            # get max PP for each site and store along with cutoff
            # (for each site, identify the cutoff with strongest signal and store to a new dataframe)

            ## END SITE-DETECTION LOOP ##

            ## BEGIN SIMULATION ROUTINES ##

            if args.sim:
                logger.info("Applying simulation-based Type I Error control...")

                # iterate over sites in the above condensed DF

                # if it is above a specified threshold (say 0.5 or 0.75)

                # get the scenario for the cutoff

                # get the profile pair for the site

                # check to see if that simulation has already been run and stored

                # if not, run it and store it

                # retrieve the requested confidence thresholds and store them in master DF

                # later, these thresholds can be mapped to Manhattan plot


                # convert simulation params to a dict
                simParams = vars(args_sim)

                # trawl the sim directory for suitable existing error curves
                cutoffsOrd, scenariosOrd, errorCurvesOrd = trawlErrorCurves(args.sim, tree_newick, scenarios, simParams)

                # generate the error curves that are missing
                for i, scenario in enumerate(scenariosOrd):
                    if not scenario:
                        scenario = newErrorCurves(tree_newick, scenario, simParams, args.sim)


        ### Collate, visualize results

        if args.heatmap or args.master_table or args.manhattan:
            # get master dataframe of results for all cutoffs
            metaDf = consolidatePCOCOutput(scenarios, scenDir, args.aa_align, logger)
            if args.key_seq:  # map the PP scores to key sequence positions
                #alignment = AlignIO.read(args.aa_align, format="fasta")
                # get key columns of the alignment
                keyCols = keyAlignmentColumns(alignment, args.key_seq)
                # filter the 2-column data table to contain only these sites
                metaDf = metaDf[metaDf["Sites"].isin(keyCols)]
                # reindex the sites
                metaDf["Sites"] = [keyCols.index(i) + 1 for i in metaDf["Sites"]]  # note index shift!
                heatmapDf = metaDf.pivot(index="cutoff", columns="Sites", values="PCOC")

        if args.heatmap:
            #print >> sys.stderr, "# MESSAGE: Drawing heatmap:"
            # Make a heatmap figure. Graphics parameters are all set in `pcoc_cont_heatmap.py`
            heatmapPath = args.output + '/' + args.heatmap
            rainbow = True  # Is the trait in question a visible wavelength to be plotted chromatically?
            heatmapper.heatMapDF(heatmapDf, heatmapPath, rainbow=rainbow)
            # tell user where it was put
            #print >> sys.stderr, heatmapPath
            logger.info("Drawing heatmap: {}".format(heatmapPath))

        if args.manhattan:

            ### Sum columns of non-exclusive probabilities
            def sumProbs(probArray):
                # the array casting is important, because numpy uses the operators element-wise
                sums = np.array([0.] * probArray.shape[1])
                for row in probArray:
                    sums = sums + np.array(row) - (sums * np.array(row))
                return sums

            #print >> sys.stderr, "# MESSAGE: Drawing Manhattan plot:"
            manhattanSeries = sumProbs(np.array(heatmapDf))  # convert to array
            manhattanPath = args.output + '/' + args.manhattan
            manhattanPlot(manhattanSeries, manhattanPath, args.key_seq)
            # tell user where it was put
            #print >> sys.stderr, manhattanPath
            logger.info("Drawing Manhattan plot: {}".format(manhattanPath))

        if args.master_table:
            #print >> sys.stderr, "# MESSAGE: Saving master table of results:"
            # Print master data table. Note that it is transpose of the table sent to heatMapDF()
            metaDf.pivot(columns="cutoff", index="Sites", values="PCOC").to_csv(args.master_table, sep='\t')
            # tell user where it was put
            #print >> sys.stderr, args.master_table
            logger.info("Saving master table of results: {}".format(args.master_table))

# try to find and load appropriate simulation results from specified directory
def trawlErrorCurves(simDir, tree_newick, scenarios, simParams):

    # set up parallel lists of cutoffs, scenarios and error curves
    # because dicts are not hashable!
    cutoffsOrd = sorted(scenarios.keys())
    scenariosOrd = [scenarios(cutoff) for cutoff in cutoffsOrd]
    errorCurvesOrd = None * len(scenariosOrd)

    # search directory for existing error curves that match and load them
    # paw thru all the files in the sim directory
    for file in os.listdir(simDir):
        curveParser = configparser.RawConfigParser()
        curveParser.read(file)
        # get pickled tree
        pickledTree = curveParser.get("sim_input", "tree")
        # get pickled scenario
        pickledScenario = literal_eval(curveParser.get("sim_input", "scenario"))
        # get pickled sim parameters
        pickledParams = dict(curveParser.items("sim_params"))

        # if parameters of the saved sim match the current user's requirements
        if ((pickledTree == tree_newick) and
                (pickledScenario in scenariosOrd) and
                (pickledParams == simParams)):
            # load it up!
            errorCurvesOrd[scenariosOrd.index(pickledScenario)] = literal_eval(
                curveParser.get("sim_results", "error_curve"))
            # and tell the user
            logger.info(
                "Error curve for cutoff {} loaded from {}.".format(cutoffsOrd[scenariosOrd.index(pickledScenario)],
                                                                   file))
        else:
            logger.debug("File {} did not match.".format(file))

    return cutoffsOrd, scenariosOrd, errorCurvesOrd

# generate a Type I Error curve for the passed tree, scenario and dict of params.
# Return the curve as dict and save a file in the specified directory.
def newErrorCurves(g_tree, simParams, simDir,
                  aProfile=None, cProfile=None, gamma=False, invGamma=False, date=None, dateTime=None,
                   falsePos=True, falseNeg=True):

    if not date:
        date = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    if not dateTime:
        dateTime = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

    # init config parser to which curve and metadata will be stored
    curveParser = configparser.RawConfigParser()
    for sname in ["sim_input", "sim_params", "sim_results"]:
        curveParser.add_section(sname)
    # store tree, .init_tree_fn should retrieve the Newick
    curveParser.set("sim_input", "tree", g_tree.init_tree_fn)
    # store scenario
    curveParser.set("sim_input", "scenario", g_tree.manual_mode_nodes)
    # store params
    for key, value in simParams:
        curveParser.set("sim_params", key, value)

    CATX_sim = simParams["CATX_sim"]
    # profiles should be passed from detection, if possible
    # but if profiles are not specified
    if not aProfile:
        aProfile = random.randrange(1, CATX_sim)
    if not cProfile:
        cProfile = aProfile
        # to ensure selection of two different profiles
        while  cProfile == aProfile:
            cProfile = random.randrange(1, CATX_sim)

    # if a Type I error curve is requested
    if falsePos:

        # generate negative alignment (all sites ancestral)
        nameAA = "%s_A%d_C%d" % (dateTime, aProfile, aProfile)
        det.bpp_lib.make_simul(nameAA, aProfile, aProfile, g_tree,
                               outputInternalSequences="no",
                               number_of_sites=simParams["n_sites"],
                               nbCAT=CATX_sim)
        fasta_file_AA = "%s/%s.fa" % (g_tree.repseq, nameAA)

        # get temp directories for pcoc_det
        tempDirs = det.get_temp_dirs(fasta_file_AA, g_tree.repseq)
        # load alignment
        ali = AlignIO.read(fasta_file_AA)

        # run site detection on sim alignment
        neg_g_tree, neg_df = det.mk_detect(g_tree.manual_mode_nodes, g_tree.init_tree_fn, ali,
                                           *tempDirs,
                                           NbCat_Est=CATX_sim,
                                           date=date,
                                           max_gap_allowed=0,  # no gaps, it's a simulated alignment!
                                           gamma=gamma,
                                           inv_gamma=invGamma
                                           )

        # generate error curve from detection df
        errorCurve = genErrorCurve(neg_df["PCOC"], errType=0b0)

        # store error curve
        curveParser.set("sim_output", "false_pos_curve", errorCurve)


    # if a Type II error curve is requested
    if falseNeg:

        # generate simulated sequences
        # positive alignment (all sites convergent)
        nameAC = "%s_A%d_C%d" % (dateTime, aProfile, cProfile)
        det.bpp_lib.make_simul(nameAC, aProfile, cProfile, g_tree,
                               outputInternalSequences = "no",
                               number_of_sites = simParams["n_sites"],
                               nbCAT = CATX_sim)
        fasta_file_AC = "%s/%s.fa" % (g_tree.repseq, nameAC)

        # get temp directories for pcoc_det
        tempDirs = det.get_temp_dirs(fasta_file_AC, g_tree.repseq)
        # load alignment
        ali = AlignIO.read(fasta_file_AC)

        # run site detection on sim alignment
        pos_g_tree, pos_df = det.mk_detect(g_tree.manual_mode_nodes, g_tree.init_tree_fn, ali,
                                       *tempDirs,
                                       NbCat_Est=CATX_sim,
                                       date=date,
                                       max_gap_allowed=0, # no gaps, it's a simulated alignment!
                                       gamma=gamma,
                                       inv_gamma=invGamma
                                       )

        # generate error curve from detection df
        errorCurve = genErrorCurve(neg_df["PCOC"], errType=0b1)

        # store error curve
        curveParser.set("sim_output", "false_neg_curve", errorCurve)


    # save error curve file
    with open(simDir + "/error_curve_" + dateTime + ".cfg") as outfile:
        curveParser.write(outfile)

    # and return error curve to caller
    return errorCurve


# take a series of PPs from a simulated alignment, and an error-type argument for whether the sites
# _should_ all call high or low (0b0 -> type I error, should all call low; 0b1 -> type II error, should all call high)
# return error rate:PP as key:value pairs
def genErrorCurve(ppSeries, errType):
    # coerce passed PPs to list and order it
    ppSeries = sorted(list(ppSeries))

    errorCurve = {}

    if errType == 0b0:
        for pp in ppSeries:
            # get the fraction of sites _above or equal to_ that PP value
            errorCurve[pp] = (100 - stat.percentileofscore(ppSeries, pp, kind="strict")) / 100
    elif errType == 0b1:
        for pp in ppSeries:
            # get the fraction of sites _below or equal to_ that PP value
            errorCurve[pp] = stat.percentileofscore(ppSeries, pp, kind="weak") / 100

    return errorCurve


### Number tree nodes for consistent reference
def init_tree(nf):
    t = Tree(nf)

    for i, n in enumerate(t.traverse("postorder")):
        n.add_features(ND = i)

    return t


### Perform BM continuous-trait ancestral state reconstruction. Return list of trait values indexed on node #.
def ancR(tree, tipTraits, logger):

    # this command does something funky to the branch lengths!
    #tree.convert_to_ultrametric(tree_length=1)

    # Check that the trait list has values for all leaves
    dataError = False
    for leaf in tree.get_leaves():
        if leaf.name not in tipTraits.keys():
            #print >> sys.stderr, "No trait data for {}!".format(leaf.name)
            logger.error("No trait data for {}!".format(leaf.name))
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
def binBy(nodeTraits, logger, fixNumBins = 0, fixBinWidth = 0):
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
        logger.info("Use of {} equal-width bins yielded {} unique trait values".format(fixNumBins, len(activeBins)))

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
        logger.info("Use of bins {} trait units wide yielded {} unique trait values".format(fixBinWidth, len(activeBins)))

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
            logger.debug("# {} <- {}".format(avg, init))
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
def uniqueScenarios(binDF, logger):

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
        logger.info("{} scenarios were consolidated into {} unique scenarios".format(numInitCols, numUniqueCols))

        for avg, init in zip(avgCutoffs, toAverage):
            #print >> sys.stderr, "# {} <- {}".format(avg, init)
            logger.debug("{} <- {}".format(avg, init))

        # construct consolidated dataframe
        binDF = pd.DataFrame(columns = avgCutoffs, data = zip(*uniqueCols))
        # sort the new averaged columns
        binDF.sort_index(axis = 1, inplace = True)

    else:

        #print >> sys.stderr, "# MESSAGE: Scenarios for all cutoffs are unique"
        logger.info("Scenarios for all cutoffs are unique.")

    return binDF


### Remove and report on scenarios in which the root is called convergent
def convergentRootFilter(binDF, logger):

    # store the original columns
    origCutoffs = binDF.columns
    # get the columns whose last (root) element is True
    convergentCutoffs = [colName for colName in origCutoffs if binDF[colName].iloc[-1]]
    # drop those columns
    binDF.drop(labels = convergentCutoffs, axis = 1, inplace = True)

    # notify user of outcome of root filter
    #print >> sys.stderr, "# MESSAGE: {}/{} scenarios eliminated due to convergent root:".format(len(convergentCutoffs), len(origCutoffs))
    logger.info("{}/{} scenarios eliminated due to convergent root: {}".format(len(convergentCutoffs), len(origCutoffs), convergentCutoffs))
    #print >> sys.stderr, "{}".format(convergentCutoffs)

    if len(convergentCutoffs) >= len(origCutoffs) / 2:
        #print >> sys.stderr, "# WARNING: Consider inverting your trait!"
        logger.warning("Consider inverting your trait!")

    return binDF

### Use tree and binary trait matrix to generate a list of 'convergent scenario' dicts usable by mk_detect()
def scenarioDicts(binDF, tree):

    # scenario strings keyed on cutoff value
    scenarios = dict()

    for col in binDF.columns:
        manual_mode_nodes = {"T": [], "C": []}
        #scenario = str()
        for node, state in enumerate(binDF[col]):
            # if convergent state is true
            if state:
                # if not the root node
                if (node < binDF.shape[0] - 1):
                    parentND = tree.search_nodes(ND=node)[0].up.ND
                    if binDF[col][parentND]:
                        manual_mode_nodes["C"].append(node)
                        #scenario = ',' + str(node) + scenario
                    else:
                        manual_mode_nodes["T"].append(node)
                        #scenario = '/' + str(node) + scenario
                else:
                    manual_mode_nodes["T"].append(node)
                    #scenario = '/' + str(node) + scenario

            scenarios[col] = manual_mode_nodes

    return scenarios


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

    return scenarios


### Gather up *all* the output files in ScenDir and put them into a master dataframe
def consolidatePCOCOutput(scenarios, scenDir, aa_align, logger):

    metaDf = pd.DataFrame()
    propTable = os.path.splitext(os.path.basename(aa_align))[:-1][
                    0] + ".results.tsv"  # get the output file from the latest run

    for cutoff in sorted(scenarios.keys()):

        try:
            subDir = scenDir + '/' + str(cutoff).replace('.','_')
            latestRun = sorted(os.listdir(subDir))[-1]
            ppFile = subDir + "/" + latestRun + "/" + propTable
            cutoffPPs = pd.read_table(ppFile)
            cutoffPPs["cutoff"] = cutoff
            metaDf = metaDf.append(cutoffPPs)
        except:
            #print >> sys.stdout, "# WARNING: Cutoff " + str(cutoff) + " not loaded"
            logger.warning("Cutoff {} not loaded".format(str(cutoff)))

    return metaDf

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

### draw a Manhattan plot of the summed PPs
def manhattanPlot(df, outPath, keyID=None):
    # at present this is meant to be run once
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    x = [i+1 for i in range(len(df))]
    y = np.array(df)
    #y = np.array([-math.log10(p) for p in y]) # log-transform
    xlimits = (min(x)-0.5, max(x)+0.5)
    plt.xlim(xlimits)
    plt.xticks(x, x, rotation=90, fontsize=7)
    plt.yticks(fontsize=20)

    if keyID:
            ax.set_xlabel("position in " + keyID, fontsize=20)
    else:
            ax.set_xlabel("alignment column")
    ax.set_ylabel("total PCOC PP", fontsize=20)

    #pThresholds = (0, 0.8, 0.9, 0.95, 1)
    #colors = ["black", "blue", "red", "violet"]
    pThresholds = (0, 0.8, 1)
    colors = ["white", (0,1,0.05)] # for black bkgd
    #colors = ["black", "blue"]  # for white bkgd
    plt.hlines(pThresholds[:-1], xlimits[0], xlimits[1], colors=colors)

    # for lin
    #masks = [[pThresholds[i] >= el > pThresholds[i + 1] for el in y] for i in range(len(pThresholds) - 1)]
    # for log
    masks = np.array([[pThresholds[i] <= el < pThresholds[i + 1] for el in y] for i in range(len(pThresholds) - 1)])

    for i, mask in enumerate(masks):
            plt.bar(x*mask, y*mask, color=colors[i], width=1)
    
    fig.set_size_inches(24, 8)
    #plt.show()
    plt.savefig(outPath)

    # no return

### the heatmap routine is imported from a separate script!

if __name__ == "__main__":
    main()