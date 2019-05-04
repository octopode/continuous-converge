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
import scipy.stats as stat
from ast import literal_eval
from ete3 import Tree#, NodeStyle, TreeStyle, TextFace, CircleFace
from shutil import rmtree
from time import sleep
from Bio import AlignIO
from collections import Counter
from pprint import pprint


# custom modules
import pcoc_treeviz as treeviz
import pcoc_plots as pviz
#from wl2rgb import wavelength_to_rgb

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
    #Options.add_argument('-i', '--invert_trait', action="store_true",
    #                     help="Invert the binary trait, i.e. assert that low trait value is the convergent state")
    Options.add_argument('-lu', '--taxa_lookup', type=str, default=None, help="name of lookup table for pretty taxon names")
    Options.add_argument('-nb', '--num_bins', type=int, default=0,
                         help="Average continuous trait into n equal-width bins, 0 = no binning")
    Options.add_argument('-bw', '--bin_width', type=float, default=0,
                         help="Average continuous trait into bins of width w, 0 = no binning")
    Options.add_argument('-m', "--min_events", type=int, default=3, help="Minimum number of transition (=convergent events). (default: 3)")
    Options.add_argument('-r', '--round_decimal', type=int, default=0,
                         help="Round traits to n decimals in filenames and figures")
    Options.add_argument('-tt', '--test_trees', action="store_true",
                         help="Draw test trees to evaluate the discretization scheme")
    Options.add_argument('-i', '--instances', type=int, default=1, help="Max number of pcoc_det instances to run concurrently (default: 1)")
    Options.add_argument('--det', action="store_true", help="Set to actually run pcoc_det.py")
    Options.add_argument('--max_gap_allowed', type=int,
                            help="max gap allowed to take into account a site (in %%), must be between 0 and 100 (default:5%%)",
                            default=5)
    Options.add_argument('--sim', action="store_true", help="Set to run post hoc simulation")
    Options.add_argument('--sim_pp_thres', type=float, default=0.5, help="Set PCOC PP threshold for running a post hoc simulation, 0.0001 = all sites")
    Options.add_argument('--sim_alpha_vals', type=float, nargs='*', default=[0.10, 0.05, 0.01], help="alpha values to test in post hoc simulations")
    Options.add_argument('--sim_beta_vals', type=float, nargs='*', default=[0.6, 0.8, 0.9], help="beta values to test in post hoc simulations")
    Options.add_argument('--sim_no_cleanup', action="store_true", help="Set to retain all the sequence sim and error estimation siles")
    Options.add_argument('-s', '--call_stationary', type=int, choices=(None,0,-1), default=None,
                         help="Override significance level of sites with stationary pairs of ML profiles. None -> bootstrap, 0 -> call indeterminate, -1 -> call negative")
    Options.add_argument('-k', '--key_seq', type=str, nargs='*', default=None, help="Names of key sequences on which to index the output columns")
    Options.add_argument('-fig', '--figure', type=int, default=None, help="Figure elements to output. 0 -> Manhattan with alignment; 1 -> Manhattan; 2 -> Shaded alignment")
    Options.add_argument('-pp', '--pp_thres', type=list, default=[0.8, 0.9, 0.95], help="PP thresholds to highlight in output")
    #Options.add_argument('-m', '--master_table', type=str, help="Save collated master data table at...")
    #Options.add_argument('-hm', '--heatmap', type=str,
    #                     help="Render heatmap from the latest available set of data and save it here")
    #Options.add_argument('-mp', '--manhattan', type=str, help="Save Manhattan plot at...")
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
                                 default=10)
    detOptions.add_argument('--gamma', action="store_true",
                                 help="Use rate_distribution=Gamma(n=4) instead of Constant()",
                                 default=False)
    detOptions.add_argument('--inv_gamma', action="store_true",
                                 help="Use rate_distribution=Gamma(n=4) instead of Constant()",
                                 default=False)
    #detOptions.add_argument('--max_gap_allowed', type=int,
    #                             help="max gap allowed to take into account a site (in %%), must be between 0 and 100 (default:5%%)",
    #                             default=5)
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

    #simOptions.add_argument('-nb_sampled_couple', type=int, metavar="INT",
    #                         help="For each convergent scenario, number of simulated alignment with different sampled couple of profiles (Ancestral/Convergent). (default: 1)",
    #                         default=1)
    simOptions.add_argument('-n_sites', type=int, metavar="INT",
                             help="Number of simulated sites per alignment. (default: 1000)",
                             default=100)
    simOptions.add_argument('-CATX_sim', type=int, choices=[10, 60],
                             help="Profile categories to simulate data (10->C10 or 60->C60). (default: set by det)",
                             default=detArgs.CATX_est)
    simOptions.add_argument('--min_dist_CAT', type=float, metavar="FLOAT",
                             help="Minimum distance between Ancestral and Convergent profiles to simulate the alignment (default: no limits)",
                             default=0)
    #simOptions.add_argument('-c_min', type=int, metavar="INT",
    #                           help="Minimum number of transition (=convergent events). (default: 2)",
    #                           default=2)
    #simOptions.add_argument('--plot_ali', action="store_true",
    #                         help="For each couple of profiles, plot a summary of the convergent scenario containing the tree and the alignment.",
    #                         default=False)
    #simOptions.add_argument('--get_likelihood_summaries', action="store_true",
    #                         help="For each couple of profiles, write a summary of the likelihoods per site.",
    #                         default=False)
    #simOptions.add_argument('--no_clean_seqs', action="store_true",
    #                         help="Do not cleanup the sequences after the run.",
    #                         default=False)
    simOptions.add_argument('-flg', type=float, metavar="FLOAT",
                               help="For each input tree, branch length multiplicator. (default: no modification)",
                               default=1)
    simOptions.add_argument('-bl_new', metavar="FLOAT", type=float,
                               help="For each input tree, replace all branch lengths by [FLOAT]. (default: no modification)",
                               default=-1)
    simOptions.add_argument('--ali_noise', action="store_true",
                               help="Add noisy events in the convergent scenario.",
                               default=False)
    simOptions.add_argument('--bl_noise', action="store_true",
                               help="Add noise in the branch lengths of of tree for the detection process.",
                               default=False)
    simOptions.add_argument('--ev_noise', choices=["+1", "-1", "=1"],
                               help="Add noise in the event placing for the detection process. +1 to add an event,  -1 to remove one, =1 to change one. -c must be fix.",
                               default=False)
    simOptions.add_argument('--root_noise', choices=["ll", "lr", "rr", "rl"],
                               help="Add noise in the root placing for the detection process. ll to move the root 2 nodes left, lr to move the root 1 node left and 1 node right, ...",
                               default=False)

    # parse the above and leave remaining argv for det
    simArgs, detArgv = simParser.parse_known_args(passthruArgv)

    return contArgs, detArgv, simArgs, simArgv


def main(contArgs, detArgv, simArgs, simArgv):
    # internal switches
    reyMethod = "PCOC"
    #reyMethod = "PC"

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
    tree_newick = tree.write()
    # not mucking with additive trees yet; ultrametrize the tree and normalize to length 1
    # 20190313: convert_to_ultrametric() is suspect. Relaxed clock with phytools instead.
    #tree.convert_to_ultrametric(tree_length=1)

    ### Generate convergent scenarios

    # Load in dict of traits keyed on species. Note hardcoding of 'sp' colName!
    tipTraits = pd.read_table(contArgs.cont_trait_table).set_index('sp').transpose().to_dict(orient='r')[0]
    # Do BM trait reconstruction
    nodeTraits = ancR(tree, tipTraits)
    # get cutoffs using specified binning
    cutoffs = binBy(nodeTraits, fixNumBins=contArgs.num_bins, fixBinWidth=contArgs.bin_width)[0]
    # Generate discrete trait matrix for all cutoff values
    #binDF = discretize(nodeTraits, cutoffs, contArgs.invert_trait)
    binDF = discretize(nodeTraits, cutoffs)
    # consolidate redundant scenarios
    # binDF = uniqueColumnFilter(binDF)[0] # element [1] is the chaff
    binDF = uniqueScenarios(binDF)
    # eliminate convergent root scenarios
    # binDF = convergentRootFilter(binDF)[0] # element [1] is the chaff
    binDF = convergentRootFilter(binDF)
    # convert binary dataframe into scenario strings
    scenarios = scenarioStrings(binDF, tree)
    # remove scenarios with too few convergent events
    scenarios = minTransitionsFilter(scenarios, contArgs.min_events)

    scenDir = os.path.abspath(contArgs.output + "/Scenarios")
    simDir = os.path.abspath(contArgs.output + "/Simulations")

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

        # get necessary number of decimals for file non-replacement
        # if specified decimal places is 0 (default; not overridden)
        if not contArgs.round_decimal:
            # get the smallest gap between cutoffs
            minGap = findMinDiff(scenarios.keys())
            decimal = -int(np.log10(minGap)) + 1
        else:
            decimal = contArgs.round_decimal

        def grayscaleMapper(x, top=max(nodeTraits), bot=min(nodeTraits), scaleMax=255):
            x = (float(x) - bot) / (top - bot)
            # low vals are white
            val = scaleMax*(1-x)
            return (val, val, val)

        # draw trait-colored tree
        treeviz.traitTree(tree, nodeTraits, grayscaleMapper, testDir, floatSwitch=decimal)
        #treeviz.traitTree(tree, nodeTraits, wavelength_to_rgb, testDir)
        # treeviz.traitTreeMinimal(tree, nodeTraits, wavelength_to_rgb, contArgs.output)  # draw trait-colored tree
        # draw boring test trees that indicate convergent scenario with red/orange
        #treeviz.testTrees(tree, scenarios, testDir, contArgs.tree, decimal)
        #treeviz.testTrees(scenarios, testDir, contArgs.tree, floatSwitch=decimal) # working 20190313
        # treeviz.testTreesMinimal(tree, scenarios, nodeTraits, wavelength_to_rgb, contArgs.tree, contArgs.output) # draw test trees with color-coded nodes
        treeviz.traitTestTrees(scenarios, nodeTraits, grayscaleMapper, testDir, contArgs.tree, floatSwitch=decimal)

    ### run site detection

    if contArgs.det:

        # make sure there are actually acceptable scenarios remaining
        if not len(scenarios):
            logger.error("There are no acceptable scenarios left to run! Adjust your criteria.")
            sys.exit(1)

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
                   "--max_gap_allowed", "100",  # always consider _all_ sites
                   "--no_cleanup"] + detArgv # add additional site-detect args from command line
        # error will be thrown if any of these are specified redundantly at the command line

        # DETECT LOOP #
        # 20190315 - with subprocessing now just compiles the commands
        logger.info("Detecting convergent sites...")

        detArgvsDynamic = [None] * len(scenarios)
        for num, cutoff in enumerate(sorted(scenarios.keys())):
            scenarioString = scenarios[cutoff]
            logger.info("convergent scenario {}/{}: cutoff {}: {}".format(num + 1, len(scenarios), cutoff, scenarioString))
            detArgvDynamic = detArgvStatic + [
                "-o", scenDir + '/' + str(cutoff).replace('.', '_'),
                "-m", str(scenarios[cutoff])]

            detArgvsDynamic[num] = detArgvDynamic

            # run pcoc.det.py
            #subprocess.call(detArgvDynamic)

        # now all the commands are ready
        subprocs = [False] * len(scenarios)
        for num, detArgvDynamic in enumerate(detArgvsDynamic):
            #while subprocs[num] is not None and subprocs[num] is False: # code for not started
            while subprocs[num] is False:
                # if the number of running processes is less than max
                if sum([safePoll(i) is None for i in subprocs]) < contArgs.instances:
                    # start up another one
                    subprocs[num] = subprocess.Popen(detArgvDynamic)
                    # this ends the while loop and moves up to the for loop to stage the next call

        # wait for all the subprocs to finish
        for subproc in subprocs:
            if isinstance(subproc, subprocess.Popen):
                subproc.wait()

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

    # get the PCOC PPs for all cutoffs
    metaDf = consolidatePCOCOutput(scenarios, scenDir, contArgs.aa_align)

    # recast the molten/long-format metaDF to wide format
    metaDfWide = metaDf.pivot(columns="Cutoff", index="Sites", values=reyMethod)

    # get the max PPs and corresponding cutoffs
    plotDf = maxPPScenarios(metaDf, methodPP = reyMethod)

    numSites = plotDf.shape[0]

    # contArgs.sim is a float value specifying how strong the a priori PP signal needs to be at a site to warrant
    # post hoc simulation.
    # this conditional block adds columns to plotDf
    if contArgs.sim:
        # make a dir for pickled sim results
        try:
            os.mkdir(simDir)
        except:
            pass
        '''
        # init post hoc simulation DF with sites index
        simDf = pd.DataFrame(index = metaDfWide.index)

        # get the highest PP for each site and the cutoff for that PP
        simDf["PP_Max"] = [max(row[1].tolist()) for row in metaDfWide.iterrows()]

        # get the cutoff for that PP
        simDf["Cutoff"] = [metaDfWide.columns[row[1].tolist().index(simDf["PP_Max"][rowNum+1])] for rowNum, row in enumerate(metaDfWide.iterrows())]
        # below should be the same command, using .iloc() for speed
        #simDf["Cutoff"] = [metaDfWide.columns[row[1].tolist().index(simDf.iloc(rowNum + 1, "PP_Max"))] for rowNum, row in enumerate(metaDfWide.iterrows())]
        '''
        # lnL matrix lookup is not by site, but by cutoff
        # so lnL matrices for all top cutoffs are loaded into memory
        # if this gets too heavy I can save/read the matrices to disk
        logger.info("Recovering ML CATegories for sites...")
        # make a set of just the relevant cutoffs
        maxPPcutoffsSet = set(plotDf.loc[plotDf["PP_Max"] >= contArgs.sim_pp_thres]["Cutoff"])
        #maxPPcutoffsSet = set(plotDf["Cutoff"])
        #TEST
        print >> sys.stderr, maxPPcutoffsSet
        # load 3D matrices for those cutoffs
        # this is a slow step, probably from all the file handling
        # its complexity should scale with number of cutoffs, not sites.
        lnLMatricesPCOC = {cutoff: consolidateBPPOutput(cutoff, scenDir, aa_align=contArgs.aa_align, withOneChange=True) for cutoff in maxPPcutoffsSet}
        lnLMatricesPC = {cutoff: consolidateBPPOutput(cutoff, scenDir, aa_align=contArgs.aa_align, withOneChange=False) for cutoff in maxPPcutoffsSet}

        # iterate over the sites and assign profiles and model lnLs to each
        '''
        for site in range(simDf.shape[0]):
            # if the site made the cut for simulation
            if simDf["PP_Max"].tolist()[site] >= contArgs.sim_pp_thres:
                # get the profiles and the model lnL and tack them on the row
                simDf.loc[site + 1, "CAT_Anc"], simDf.loc[site + 1, "CAT_Con"], simDf.loc[site + 1, "lnL_PCOC"] = getMLCATProfiles(site, lnLMatricesPCOC[simDf["Cutoff"].tolist()[site]])
                #simDf.loc[site + 1, "CAT_Anc_PC"], simDf.loc[site + 1, "CAT_Con_PC"], simDf.loc[site + 1, "lnL_PC"] = getMLCATProfiles(site, lnLMatricesPC[cutoff])
        '''
        for site in range(plotDf.shape[0]):
            # if the site made the cut for simulation
            if plotDf["PP_Max"].tolist()[site] >= contArgs.sim_pp_thres:
                # get the profiles and the model lnL and tack them on the row
                if reyMethod == "PCOC":
                    plotDf.loc[site + 1, "CAT_Anc"], plotDf.loc[site + 1, "CAT_Con"], plotDf.loc[site + 1, "lnL_PCOC"] = getMLCATProfiles(site, lnLMatricesPCOC[plotDf["Cutoff"].tolist()[site]])
                elif reyMethod == "PC":
                    plotDf.loc[site + 1, "CAT_Anc"], plotDf.loc[site + 1, "CAT_Con"], plotDf.loc[site + 1, "lnL_PC"] = getMLCATProfiles(site, lnLMatricesPC[plotDf["Cutoff"].tolist()[site]])
                    #plotDf.loc[site + 1, "CAT_Anc_PC"], plotDf.loc[site + 1, "CAT_Con_PC"], plotDf.loc[site + 1, "lnL_PC"] = getMLCATProfiles(site, lnLMatricesPC[cutoff])

        # display profile numbers as ints in the table
        # but this crashes when there are NaNs :/
        #simDf["CAT_Anc"] = simDf["CAT_Anc"].astype(int)
        #simDf["CAT_Con"] = simDf["CAT_Con"].astype(int)

        #postHocFilename = os.path.splitext(contArgs.master_table)[0] + "_posthoc.tsv"
        #simDf.to_csv(postHocFilename, sep='\t')

        #TEST
        #print plotDf

        # make a list of the different error curves required
        # this amounts to the unique rows of simDf, less the lnL column
        #uniqueSims = simDf.dropna()[["Cutoff", "CAT_Anc", "CAT_Con"]].drop_duplicates().reset_index(drop = True)
        uniqueSims = plotDf.dropna()[["Cutoff", "CAT_Anc", "CAT_Con"]].drop_duplicates().reset_index(drop=True)
        # Not sure whether this is necessary: convert the profiles to ints for merging (the NaNs are gone!)
        uniqueSims["CAT_Anc"] = uniqueSims["CAT_Anc"].astype(int)
        uniqueSims["CAT_Con"] = uniqueSims["CAT_Con"].astype(int)
        # add a column of scenario strings
        uniqueSims["Scenario"] = [scenarios[cutoff] for cutoff in uniqueSims["Cutoff"].tolist()]

        #print "unique required sims"
        #print uniqueSims

        empty = contArgs.output + "/empty"

        pickleDir = os.path.join(simDir, "Pickles")

        # make pickle dir
        if not os.path.isdir(pickleDir):
            os.mkdir(pickleDir)

        # load all the available error curves into a DataFrame: metadata mapped to list of confidence thresholds
        availSims = loadSavedSims(pickleDir, tree_newick, simArgs, alpha = contArgs.sim_alpha_vals, beta = contArgs.sim_beta_vals)
        #print "available sims"
        #print availSims

        # left-join those to the table of unique required sims
        uniqueSims = pd.merge(uniqueSims, availSims, how='left', on=["Scenario", "CAT_Anc", "CAT_Con"])
        # and if there are multiple pickles appropriate for each scenario, use only the freshest ones
        uniqueSims.drop_duplicates(subset = ["Scenario", "CAT_Anc", "CAT_Con"], keep = "last", inplace = True)
        #print "sims recovered from picklejar"
        #print uniqueSims

        # get rows for sims that still need to be run
        newSims = uniqueSims.loc[uniqueSims["PP_Threshold_TypeI"].isnull()]
        # run them and update the table with confidence thresholds
        #print "missing sims"
        newSims = runSiteSims(newSims, contArgs, simArgv, simDir, pickleDir, simArgs, reyMethod=reyMethod, alpha = contArgs.sim_alpha_vals, beta = contArgs.sim_beta_vals)

        # replace NaNs in threshold cols with empty list
        for row in uniqueSims.loc[uniqueSims.PP_Threshold_TypeI.isnull(), 'PP_Threshold_TypeI'].index:
            uniqueSims.at[row, 'PP_Threshold_TypeI'] = []
        for row in uniqueSims.loc[uniqueSims.PP_Threshold_TypeII.isnull(), 'PP_Threshold_TypeII'].index:
            uniqueSims.at[row, 'PP_Threshold_TypeII'] = []
        
        uniqueSims = mergeFillNans(uniqueSims, newSims, how='left', on=["Cutoff", "CAT_Anc", "CAT_Con", "Scenario"])

        # left-join the unique sims table to the sitewise table
        #plotDf = pd.merge(simDf, uniqueSims, how='left', on=["Cutoff", "CAT_Anc", "CAT_Con"])
        plotDf = pd.merge(plotDf, uniqueSims, how='left', on=["Cutoff", "CAT_Anc", "CAT_Con"])

        # TEST
        #print plotDf

    # ascertain significance level for each site
    # for bootstrap, it goes ...-3, -2, -1, 0, 1, 2, 3... where -3, 3 are the most stringent negative and positive thres
    if contArgs.sim:
        plotDf["SigLevel"] = np.vectorize(rankPosNeg)(plotDf['PP_Max'], plotDf['PP_Threshold_TypeI'],
                                                      plotDf['PP_Threshold_TypeII'])
    else:
    # for PP, it goes 0, 1, 2, 3 etc for the specified cutoffs
        plotDf['PP_Threshold'] = [contArgs.pp_thres] * numSites
        plotDf["SigLevel"] = np.vectorize(rank)(plotDf["PP_Max"], plotDf['PP_Threshold'])

    # get gap percentage for each site in alignment
    ali = AlignIO.read(contArgs.aa_align, "fasta")
    plotDf["GapPct"] = gapPcts(ali)

    if contArgs.sim:
        # check conservation
        plotDf["ConservedPct"] = consPcts(ali)

        # designate unacceptably gappy sites as indeterminate
        logger.debug("Silencing gappy sites...")
        for i in range(plotDf.shape[0]):
            # too gappy?
            if plotDf.loc[i, "GapPct"] > contArgs.max_gap_allowed:
                plotDf.at[i, "PP_Max"] = np.nan
                plotDf.at[i, "PP_Threshold_TypeI"] = []
                plotDf.at[i, "PP_Threshold_TypeII"] = []
                plotDf.at[i, "SigLevel"] = 0
            # check that fully conserved sites call properly
            if plotDf.loc[i, "ConservedPct"] == 100 and plotDf.loc[i, "SigLevel"] > -1:
                # issue a warning if they do not
                logger.debug("Site {} is fully conserved, but calls with significance level {}.\nConsider using more bootstrap replicates."
                               .format(i+1, plotDf.loc[i, "SigLevel"]))
                plotDf.at[i, "SigLevel"] = -1

        # designate apparently stationary sites as indeterminate or negative (as specified)
        if contArgs.call_stationary is not None:
            for i in range(plotDf.shape[0]):
                if plotDf.loc[i, "CAT_Anc"] == plotDf.loc[i, "CAT_Con"]:
                    #plotDf.at[i, "PP_Max"] = np.nan
                    plotDf.at[i, "PP_Threshold_TypeI"] = []
                    plotDf.at[i, "SigLevel"] = contArgs.call_stationary
                    if contArgs.call_stationary == 0:
                        plotDf.at[i, "PP_Threshold_TypeII"] = []
                    elif contArgs.call_stationary < 0:
                        plotDf.at[i, "PP_Threshold_TypeII"] = [plotDf.loc[i, "PP_Max"]] + ([0] * (len(contArgs.sim_beta_vals) - 1))

    # get basename of the alignment to name outfiles after
    # should change this to basename of the trait file, since doing multiple traits
    aliBasename = os.path.splitext(os.path.basename(contArgs.aa_align))[0]

    # TEST
    print plotDf

    # map the PP scores to key sequence positions if desired
    # alignment reindexing is done post-analysis so user can switch key seqs on the fly without reanalyzing
    if contArgs.key_seq:
        # read in the alignment
        ali = AlignIO.read(contArgs.aa_align, format="fasta")
        # for each key species
        for sp in contArgs.key_seq:
            # get key columns of the alignment
            #keyCols = keyAlignmentColumns(ali, contArgs.key_seq)
            ## filter the data table to contain only these sites
            #plotDfMapped = plotDf[plotDf["Sites"].isin(keyCols)]
            ## reindex the sites from 1 to n
            #plotDfMapped["Sites"] = range(1, len(plotDfMapped.index) + 1)
            plotDfMapped = reindexDataframe(plotDf, ali, sp)

            # add Sites index
            plotDfMapped["Sites"] = range(1, plotDfMapped.shape[0] + 1)
            plotDfMapped.set_index("Sites", inplace=True)

            # save results table keyed on sp
            tsvPath = os.path.join(contArgs.output, aliBasename + "_key-" + sp + "_results.tsv")
            plotDfMapped.to_csv(tsvPath, sep='\t')
            # tell user where it was put
            logger.info("Table of results saved at {}".format(tsvPath))

            # if figures are requested
            if contArgs.figure is not None:
                # see if there is a lookup table for pretty taxon names
                if contArgs.taxa_lookup:
                    # Load in dict of traits keyed on species. Note hardcoding of 'sp' colName!
                    luTable = pd.read_table(contArgs.taxa_lookup, header=False, sep='\t')
                    prettySeqNames = luTable.set_index(luTable.columns[0]).transpose().to_dict(orient='r')[0]
                else:
                    prettySeqNames = None

                figPath = os.path.join(contArgs.output, aliBasename + "_key-" + sp + "_figure.pdf")
                # remap the alignment itself
                aliMapped = reindexAlignment(ali, sp)
                # save figure keyed on sp
                pviz.masterFigure(plotDfMapped, aliMapped, tipTraits, elements=contArgs.figure, prettySeqNames=prettySeqNames,
                                  outPath=figPath, alpha=contArgs.sim_alpha_vals, beta=contArgs.sim_beta_vals, xLabel="amino acid site in {}".format(sp))
                logger.info("Figure saved at {}".format(figPath))

    # results are always output for the original alignment, with gaps

    # output collated data. This is no longer subject to an a command-line arg
    # Print master data table. Note that it is transpose of the table sent to heatMapDF()
    # metaDfWide.to_csv(contArgs.master_table, sep='\t')
    tsvPath = os.path.join(contArgs.output, aliBasename + "_results.tsv")

    # add Sites index
    plotDf["Sites"] = range(1, plotDf.shape[0] + 1)
    plotDf.set_index("Sites", inplace=True)

    plotDf.to_csv(tsvPath, sep='\t')
    # tell user where it was put
    logger.info("Table of results saved at {}".format(tsvPath))

    # draw the requested kind of figure
    if contArgs.figure is not None:
        figPath = os.path.join(contArgs.output, aliBasename + "_figure.pdf")
        #ali = AlignIO.read(contArgs.aa_align, "fasta") # moved up for gap checking
        # draw master figure
        pviz.masterFigure(plotDf, ali, tipTraits, elements=contArgs.figure, outPath=figPath, alpha=contArgs.sim_alpha_vals, thresPP=contArgs.pp_thres, beta=contArgs.sim_beta_vals)
        logger.info("Figure saved at {}".format(figPath))

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

### Invert scenarios in which the root is called convergent
def convergentRootFilter(binDF):

    # store the original columns
    origCutoffs = binDF.columns
    # get the columns whose last (root) element is True
    convergentCutoffs = [colName for colName in origCutoffs if binDF[colName].iloc[-1]]
    # invert those columns
    for cutoff in convergentCutoffs:
        binDF[cutoff] = ~binDF[cutoff]
    # drop those columns
    #binDF.drop(labels = convergentCutoffs, axis = 1, inplace = True)

    # notify user of outcome of root filter
    logger.info("{}/{} scenarios inverted due to convergent root:\n{}".format(len(convergentCutoffs), len(origCutoffs), convergentCutoffs))

    return binDF

### Remove and report on scenarios in which the root is called convergent
def convergentRootFilterOld(binDF):

    # store the original columns
    origCutoffs = binDF.columns
    # get the columns whose last (root) element is True
    convergentCutoffs = [colName for colName in origCutoffs if binDF[colName].iloc[-1]]
    # drop those columns
    binDF.drop(labels = convergentCutoffs, axis = 1, inplace = True)

    # notify user of outcome of root filter
    logger.info("{}/{} scenarios eliminated due to convergent root:\n{}".format(len(convergentCutoffs), len(origCutoffs), convergentCutoffs))

    if len(convergentCutoffs) >= len(origCutoffs) / 2:
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


### Remove and report on scenarios that have less than the specified number of convergent events
def minTransitionsFilter(scenarios, minEvents):

    wheat = {key: value for key, value in scenarios.iteritems() if (value.count('/') + 1) >= minEvents}
    chaff = {key: value for key, value in scenarios.iteritems() if (value.count('/') + 1) < minEvents}

    # notify user of filter action
    logger.info("{}/{} scenarios eliminated due to fewer than {} convergent events:\n{}".format(len(chaff), len(scenarios), minEvents, sorted(chaff.keys())))

    return wheat


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

    # for each cutoff
    for cutoff in sorted(scenarios.keys()):
        found = False
        subDir = scenDir + '/' + str(cutoff).replace('.', '_')
        # paw thru the stored results for the most recent set with matching alignment name
        for runDir in sorted(os.listdir(subDir), reverse=True):
            #latestRun = sorted(os.listdir(subDir))[-1]
            # if the results file matches the alignment
            if ppFileBasename in os.listdir(os.path.join(subDir, runDir)):
                found = True
                # load it in
                #ppFile = subDir + "/" + runDir + "/" + ppFileBasename
                ppFile = os.path.join(subDir, runDir, ppFileBasename)
                logger.info("Loading PPs for cutoff {} from {}.".format(cutoff, ppFile))
                cutoffPPs = pd.read_table(ppFile)
                # store the cutoff and the scenario string
                cutoffPPs["Cutoff"] = cutoff
                cutoffPPs["Scenario"] = scenarios[cutoff]
                metaDf = metaDf.append(cutoffPPs)
                break
        # if the search is fruitless
        if not found:
            logger.warning("Posterior probabilities not loaded for cutoff {}.".format(cutoff))

    return metaDf


# take the dataframe returned by consolidatePCOCOutput()# and return a sitewise DF with columns `Sites`, `Cutoff`, `PP_max`
def maxPPScenarios(metaDf, methodPP="PCOC"):

    # recast the molten/long-format metaDF to wide format
    metaDfWide = metaDf.pivot(columns="Cutoff", index="Sites", values=methodPP)

    # init new DF with sites index
    plotDf = pd.DataFrame(index=metaDfWide.index)

    # iterrows() calls could probably be optimized
    # get the highest PP for each site and the cutoff for that PP
    plotDf["PP_Max"] = [max(row[1].tolist()) for row in metaDfWide.iterrows()]
    # get the cutoff for that PP
    plotDf["Cutoff"] = [metaDfWide.columns[row[1].tolist().index(plotDf["PP_Max"][rowNum + 1])] for rowNum, row in
                       enumerate(metaDfWide.iterrows())]

    return plotDf


### Gather up the .infos files for a particular cutoff in ScenDir and put the lnL data into a 3D array
def consolidateBPPOutput(cutoff, scenDir, aa_align, withOneChange = True):

    # search criterion for files to load from the Estimations dir
    if withOneChange:
        searchPattern = "withOneChange.infos"
    else:
        searchPattern = "noOneChange.infos"

    # get the name of the output file
    ppFileBasename = os.path.splitext(os.path.basename(aa_align))[:-1][0] + ".results.tsv"
    # name of subdirectory for applicable cutoff
    subDir = scenDir + '/' + str(cutoff).replace('.', '_')
    found = False

    # paw thru the stored results for the most recent set with matching alignment name
    for runDir in sorted(os.listdir(subDir), reverse=True):
        # run folders are named by numeric datetime, so last is most recent
        #latestRun = sorted(os.listdir(subDir))[-1]
        # if the results file matches the alignment
        if ppFileBasename in os.listdir(os.path.join(subDir, runDir)):
            found = True
            #estimsDir = subDir + "/" + runDir + "/Estimations"
            estimsDir = os.path.join(subDir, runDir, "Estimations")
            logger.info("Recovering ML CATegories for cutoff {} from {}.".format(cutoff, estimsDir))
            # get list of target filenames
            lnLfnames = sorted([fname for fname in os.listdir(estimsDir) if fname.endswith(searchPattern)])
            break

    if not found:
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
    #lnLs = np.empty((numProfilesA, numSites, numProfilesC))
    #lnLs = np.empty((numProfilesA, numProfilesC, numSites)) # seems right?

    # iteratively build the array
    for i, fname in enumerate(lnLfnames):
        #print pd.read_table(estimsDir + "/" + lnLfnames[i])["lnL"].tolist()
        # more than one non-: entry in the slicing tuple?
        #lnLs[:][profilesA[i]-1][profilesC[i]-1] = pd.read_table(estimsDir + "/" + lnLfnames[i])["lnL"].tolist()
        lnLs[:, profilesA[i] - 1, profilesC[i] - 1] = pd.read_table(estimsDir + "/" + lnLfnames[i])["lnL"].tolist()

    return lnLs

# get the ancestral and convergent CAT profiles from name of a PCOC/BPP .infos file
def infosFname2Profiles(fname):
    # split filename on '_'; take the 2nd- and 3rd-to-last elements
    return [int(cat) for cat in fname.split('_')[-3:-1]]

### Take a site (0-indexed) and a np.array of log-likelihoods for different ancestral/convergent CAT models,
### return the ancestral and convergent CATegory pair with the highest log likelihood. Also return lnL for diagnostics.
def getMLCATProfiles(site, lnLs):

    # get A, C indices of the max lnL for a given site in the matrix
    # This is the only spot where the ordering of index could switch up A, C.
    # 20190315: I checked its behavior and it is correct!
    profileA, profileC = np.unravel_index(np.argmax(lnLs[site, :, :]), (lnLs.shape[1], lnLs.shape[2]))
    # get the actual lnL
    lnL = lnLs[site, profileA, profileC]

    # return 1-indexed CAT numbers
    return profileA + 1, profileC + 1, lnL

# search directory for all pickled error curves generated using the specified .__dict__ of args
def loadSavedSims(pickleDir, tree_newick, simArgs, alpha = [0.10, 0.05, 0.01], beta = [0.6, 0.8, 0.9]):

    # paw thru all the pickles in the jar
    files = os.listdir(pickleDir)
    # init list of DF rows. Cannot know length b/c don't know how many of the sim files are valid.
    rows = list()
    for file in files:
        curveParser = configparser.RawConfigParser()
        for sname in ["sim_args", "sim_input", "sim_output"]:
            curveParser.add_section(sname)
        curveParser.readfp(open(os.path.join(pickleDir, file), 'r'))

        # check that all pickled pcoc_sim [static] args are same or better than specified
        # this should take a dict coming from parse_args(), not the argv that actually went to the pcoc_sim call
        pickleCheckList = list()
        for key, value in simArgs.__dict__.items():
            try:
                pickleCheckList.append(curveParser.getint("sim_args", key) >= value)
            except:
                # if one of the passed options is not stored in the pickle
                logger.debug("{} is not in pickle {}!\nAnalysis proceeds, but consider stashing your pickles.".format(key, file))
                pickleCheckList.append(True)
                pass

        if (all(pickleCheckList) and curveParser.get("sim_input", "tree_newick") == tree_newick):
            # init DataFrame row
            row = pd.Series()
            # load it up! (from dict)
            errorCurveTypeI = literal_eval(curveParser.get("sim_output", "error_curve_type1"))
            # get the PP values corresponding to passed levels of alpha
            row.at["PP_Threshold_TypeI"] = getConfThresholds(errorCurveTypeI, alpha)
            errorCurveTypeII = literal_eval(curveParser.get("sim_output", "error_curve_type2"))
            # get the PP values corresponding to passed levels of alpha
            row.at["PP_Threshold_TypeII"] = sorted(getConfThresholds(errorCurveTypeII, beta), reverse = True)
            # get dynamic sim metadata
            row["Scenario"] = curveParser.get("sim_input", "scenario")
            row["CAT_Anc"] = curveParser.getint("sim_input", "profile_anc")
            row["CAT_Con"] = curveParser.getint("sim_input", "profile_con")
            rows.append(row.to_frame(1).T)

    # concat the series into a DF
    if len(rows):
        simsCatalog = pd.concat(rows, ignore_index = True)
    else:
        # empty default DF
        simsCatalog = pd.DataFrame(columns=["Scenario", "CAT_Anc", "CAT_Con", "PP_Threshold_TypeI", "PP_Threshold_TypeII"])

    logger.info("Loaded {} pickled simulations from {}.".format(simsCatalog.shape[0], pickleDir))

    return simsCatalog

# run simulations as requested in the passed pd.DataFrame()
def runSiteSims(newSims, contArgs, simArgv, simDir, pickleDir, simArgs, reyMethod="PCOC", alpha = [0.10, 0.05, 0.01], beta = [0.6, 0.8, 0.9]):

    dateTime = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    tree = init_tree(contArgs.tree)
    tree_newick = tree.write()

    #["Scenario", "CAT_Anc", "CAT_Con", "PP_Threshold_TypeI"]
    # location of the modified version of pcoc_sim.py
    # it accepts manual profile specification and a single tree path
    cwd = os.path.abspath(os.path.join(os.path.realpath(__file__), '..'))
    modifiedPCOCSim = os.path.join(cwd, "pcoc_sim_manual.py")

    # set args that are the same for all sims
    simArgvStatic = [modifiedPCOCSim,
                     "-td", contArgs.tree,
                     #"-o", simDir,
                     #"--pcoc",
                     "--no_clean_seqs"] + simArgv  # add additional site-detect args from command line
    # error will be thrown if any of these are specified redundantly at the command line

    # SIMULATE LOOP #

    newSims.reset_index(inplace = True)
    numSims = newSims.shape[0]
    logger.info("Running {} new post hoc simulations...".format(numSims))

    # a list of sim-det processes only
    subprocs = [False] * numSims * 2
    for i, row in newSims.iterrows():
        # keep checking until the process is started
        while subprocs[i * 2] is False:
            # if the number of running processes is less than max
            if sum([safePoll(j) is None for j in subprocs]) < contArgs.instances:
                profileA = int(row["CAT_Anc"])
                profileC = int(row["CAT_Con"])
                scenStr = row["Scenario"]
                simIDStr = "A{}_C{}_{}".format(profileA, profileC, scenStr)
                logger.info("Simulation {}/{}: {}".format(i+1, numSims, simIDStr))

                # sanitize simulation identifier for use as a directory name
                # it might later be necessary to just call the directory "Sim_i"
                # if the dirnames get too long
                #simIDStr = simIDStr.replace('/', '_').replace(',', "-")
                simIDStr = "sim_{}_{}".format(i, dateTime)
                # make output dir
                simOutDir = os.path.join(simDir, simIDStr)
                try:
                    os.mkdir(simOutDir)
                except:
                    pass

                simArgvDynamic = simArgvStatic + [
                    "-o", simOutDir,
                    "-m", str(row["Scenario"]),
                    "-m_sampled_couple", str(profileA), str(profileC),
                    ]

                # simulate the alignment
                # use call() here because it needs to finish before the other 2 start
                subprocess.call(simArgvDynamic)

                # a clumsy way to make sure detect results don't end up in the sim directory:
                sleep(2)  # since they are named using dateTime

                # detect on both NEGATIVE, POSITIVE algts for Type I, II error
                detArgvStatic = ["pcoc_det.py", # program call
                                 "-f", "-1", # don't waste time filtering output
                                 ]

                # iterate over negative, then positive
                for type, profileC in enumerate((int(row["CAT_Anc"]), int(row["CAT_Con"]))):
                    # if we're dealing with a stationary site
                    if profileA == profileC and type == 1:
                        # symlink the pos dir to the neg dir
                        #os.symlink(os.path.join(simOutDir, "neg"), os.path.join(simOutDir, "pos"))
                        # the detection has already been done, so set job status to done and move on
                        subprocs[i * 2 + type] = True
                    else:
                        # retrieve the simulated alignment
                        simAlgtPath = os.path.join(simOutDir, sorted(os.listdir(simOutDir))[0], "Tree_1/sequences/Scenario_1",
                                                   "Scenario_1_A{}_C{}.fa".format(profileA, profileC))

                        # make _another_ subdirectory so the det jobs don't trip over each other
                        if type == 0:
                            detOutDir = os.path.join(simOutDir, "neg")
                            os.mkdir(detOutDir)
                        elif type == 1:
                            detOutDir = os.path.join(simOutDir, "pos")
                            os.mkdir(detOutDir)

                        # set up detection on the simulated alignment
                        detArgvDynamic = detArgvStatic + [
                                   "-t", contArgs.tree,
                                   "-aa", simAlgtPath,
                                   "-o", detOutDir,
                                   "-m", str(row["Scenario"]),
                                         ] + detArgv # add additional site-detect args from command line
                        # error will be thrown if any of these are specified redundantly at the command line

                        # store the Popen object for polling
                        subprocs[i * 2 + type] = subprocess.Popen(detArgvDynamic)
                        # wait 5 seconds before trying to spawn any more subprocesses
                        sleep(2)

    # wait for all the subprocs to finish
    for subproc in subprocs:
        if isinstance(subproc, subprocess.Popen):
            subproc.wait()

    for i, row in newSims.iterrows():

        # recreate from DF all this info used above:
        profileA = int(row["CAT_Anc"])
        profileC = int(row["CAT_Con"])
        scenStr = row["Scenario"]
        #logger.info("Simulation {}/{}: {}".format(i + 1, numSims, simIDStr))

        # sanitize simulation identifier for use as a directory name
        # it might later be necessary to just call the directory "Sim_i"
        # if the dirnames get too long
        # simIDStr = simIDStr.replace('/', '_').replace(',', "-")
        simIDStr = "sim_{}_{}".format(i, dateTime) # dateTime is static within this function
        # make output dir
        simOutDir = os.path.join(simDir, simIDStr)

        # retrieve FPR results table (2nd run)
        detOutDir = os.path.join(simOutDir, "neg")
        simResultsTypeIPath = os.path.join(simOutDir, detOutDir, os.listdir(detOutDir)[0], "Scenario_1_A{}_C{}.results.tsv".format(profileA, profileA))
        simResultsTypeI = pd.read_table(simResultsTypeIPath)
        # generate an error curve from it
        errorCurveTypeI = genErrorCurve(simResultsTypeI[reyMethod], 0b0)
        # if not a stationary sim
        if profileA != profileC:
            # retrieve FNR results table
            detOutDir = os.path.join(simOutDir, "pos")
            simResultsTypeIIPath = os.path.join(simOutDir, detOutDir, os.listdir(detOutDir)[0], "Scenario_1_A{}_C{}.results.tsv".format(profileA, profileC))
            simResultsTypeII = pd.read_table(simResultsTypeIIPath)
            # generate an error curve from it
            errorCurveTypeII = genErrorCurve(simResultsTypeII[reyMethod], 0b1)
        # if it is stationary, type II error curve makes no sense
        else:
            errorCurveTypeII = {1:1}
        # and pickle it in a file
        picklePath = os.path.join(pickleDir, simIDStr + ".conf")
        pickleErrorCurve(errorCurveTypeI, errorCurveTypeII, simArgs, tree_newick, (profileA, profileC), row["Scenario"], picklePath)
        # get confidence thresholds and add them to newSims DF
        #newSims.iloc[i, newSims.columns.get_loc("PP_Threshold_TypeI")] = getConfThresholds(errorCurveTypeI, alpha)
        # .at[] allows insertion of a list into DF
        newSims.at[i, "PP_Threshold_TypeI"] = getConfThresholds(errorCurveTypeI, alpha)
        newSims.at[i, "PP_Threshold_TypeII"] = sorted(getConfThresholds(errorCurveTypeII, beta), reverse = True)
        if not contArgs.sim_no_cleanup:
            # delete all sim files besides the error curve
            rmtree(simOutDir)

    # END SIMULATE LOOP #

    # return the completed sim table
    return newSims

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

def pickleErrorCurve(errorCurveTypeI, errorCurveTypeII, simArgs, tree_newick, profiles, scenario, filename):

    # init config parser to which curve and metadata will be stored
    curveParser = configparser.RawConfigParser()
    for sname in ["sim_args", "sim_input", "sim_output"]:
        curveParser.add_section(sname)

    # store sim args (including default ones)
    for key, value in simArgs.__dict__.items():
        curveParser.set("sim_args", key, value)

    # store tree Newick
    curveParser.set("sim_input", "tree_newick", tree_newick)
    # store convergent scenario
    curveParser.set("sim_input", "scenario", scenario)
    # store profiles
    curveParser.set("sim_input", "profile_anc", profiles[0])
    curveParser.set("sim_input", "profile_con", profiles[1])

    # store error curve
    curveParser.set("sim_output", "error_curve_type1", errorCurveTypeI)
    curveParser.set("sim_output", "error_curve_type2", errorCurveTypeII)

    with open(filename, 'w') as handle:
        curveParser.write(handle)

    return filename

# take an error curve (dict form) and set of alpha values
# return the PP thresholds for those alpha values
def getConfThresholds(errorCurve, errorRates):

    # sort the alpha values descending
    errorRates = sorted(errorRates, reverse=True)
    ppThresholds = [None] * len(errorRates)
    i = 0
    for key in sorted(errorCurve.keys()):
        if i < len(errorRates) and errorCurve[key] <= errorRates[i]:
            ppThresholds[i] = key
            i += 1

    return ppThresholds


def reindexDataframe(df, ali, keyID):

    # get the key sequence
    seqIDs = [record.id for record in ali]
    keySeq = ali[seqIDs.index(keyID)].seq
    gapIndices = list(findOffsets(keySeq, '-'))

    # get the DF index labels corresponding to those integer positions
    gapIndexLabels = [df.index[indexPos] for indexPos in gapIndices]
    # drop corresponding rows from DF
    # making a copy so as not to mess up DF for calling process
    outDf = df.drop(gapIndexLabels, axis=0)
    outDf.reset_index(drop=True, inplace=True)

    # reset the "fake" site index if it exists
    if "Sites" in outDf.columns:
        outDf["Sites"] = [i + 1 for i in outDf.index]

    return outDf

### take a BioAlign Alignment and return it with gaps in the keyID seq removed
def reindexAlignment(ali, keyID):
    # get the key sequence
    seqIDs = [record.id for record in ali]
    keySeq = ali[seqIDs.index(keyID)].seq
    aliLength = ali.get_alignment_length()

    aliSlices = []
    for i in range(aliLength):  # only iterate over columns that are not gaps in target seq
        if keySeq[i] != "-":  # meaning anything except gaps
            aliSlices.append(ali[:, i: i+1])

    return sum(aliSlices, ali[:,0:0])

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
            mappedCols[keyIndex] = i  # map the alignment column to the key column.
            keyIndex += 1

    return mappedCols

### Find the start of all (possibly-overlapping) instances of needle in haystack
# this generator from https://stackoverflow.com/questions/11122291/python-find-char-in-string-can-i-get-all-indexes
def findOffsets(haystack, needle):

    offs = -1
    while True:
        offs = haystack.find(needle, offs+1)
        if offs == -1:
            break
        else:
            yield offs

### Utility function that does a merge, then consolidates duplicated columns.
### It does this by choosing the first truthy value in the pair of columns for each row.
def mergeFillNans(left, right, **kwargs):

    '''
    # a helper function that makes np.nan Falsy and gives single value for truth of lists
    def safeBool(x):
        if not isinstance(x, list):
            if pd.isnull(x): return False
            else:return bool(x)
        else:
            return any(x)
    '''

    merged = pd.merge(left, right, **kwargs).astype('object')
    dupedColumns = [colname[:-2] for colname in merged.columns if "_x" in colname]

    # pick the max values along the whole pair of columns (np.nan = -inf)
    for colname in dupedColumns:
        colnamex = colname + "_x"
        colnamey = colname + "_y"
        colx = merged[colnamex].tolist()
        coly = merged[colnamey].tolist()
        colMerged = [max(comp) for comp in zip(colx, coly)]
        merged[colname] = colMerged
        merged.drop([colnamex, colnamey], axis=1, inplace=True)

    return merged

# Take Bio.Align object, return a list of the percentage of gaps in each column
def gapPcts(algt):
    return [(100 * float(str(algt[:, site]).count('-')) / float(len(algt))) for site in range(algt.get_alignment_length())]

# Take Bio.Align object, return a list of the percent conservation of each column (rate of the most common AA)
def consPcts(algt):
    return [(100 * float(Counter(str(algt[:, site])).most_common(1)[0][1]) / float(len(algt))) for site in range(algt.get_alignment_length())]

# Returns minimum difference between any pair
def findMinDiff(arr):

    n = len(arr)
    # Sort array in non-decreasing order
    arr = sorted(arr)

    # Initialize difference as infinite
    diff = 10**20

    # Find the min diff by comparing adjacent
    # pairs in sorted array
    for i in range(n-1):
        if arr[i+1] - arr[i] < diff:
            diff = arr[i+1] - arr[i]

    # Return min diff
    return float(diff)

### Poll return status of a Popen object
# if not a Popen object, return the argument unchanged
def safePoll(pop):
    if isinstance(pop, subprocess.Popen):
        return pop.poll()
    else:
        return pop

### Rank val in iterable seri
# has a built-in sort
def rank(val, seri):
    if isinstance(seri, list):
        rank = 0
        for i, thres in enumerate(sorted(seri)):
            if val >= thres:
                rank += 1
        return rank
        # should return 0 if len(seri) == 0
    else:
        return int(val >= seri)

# rank() wrapper that takes positive and negative cutoffs
# and returns a signed sig level
def rankPosNeg(val, seriPos, seriNeg):

    # check if it's positive
    sigLevelPos = rank(val, seriPos)

    # check if it's negative
    if isinstance(seriNeg, list):
        sigLevelNeg = rank(val, seriNeg) - len(seriNeg)
    # if only one negative sigLevel
    else:
        sigLevelNeg = rank(val, seriNeg) - 1

    # if PP is in the overlap zone or the gap
    if (sigLevelPos > 0 and sigLevelNeg < 0) or (sigLevelPos == sigLevelNeg == 0):
        sigLevel = 0
    else:
        # return the one that is not 0!
        sigLevel = [i for i in (sigLevelPos, sigLevelNeg) if i != 0][0]

    return sigLevel

### Number tree nodes for consistent reference
def init_tree(nf):
    t = Tree(nf)

    for i, n in enumerate(t.traverse("postorder")):
        n.add_features(ND = i)

    return t

if __name__ == "__main__":
    contArgs, detArgv, simArgs, simArgv = parse_args(sys.argv)
    main(contArgs, detArgv, simArgs, simArgv)