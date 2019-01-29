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

import os
import sys
import subprocess
import argparse
import time
import numpy as np
import pandas as pd
from ete3 import Tree, NodeStyle, TreeStyle, TextFace
from matplotlib.colors import LinearSegmentedColormap
from re import findall, search
from pprint import pprint
from Bio import AlignIO
from wl2rgb import wavelength_to_rgb
import pcoc_cont_heatmap as heatmapper
import pgls2_3 as cpgls

##########
# inputs #
##########
start_time = time.time()

### Option defining
parser = argparse.ArgumentParser(prog="pcoc_num_tree.py",
                                 description='')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

##############
requiredOptions = parser.add_argument_group('Required arguments')
requiredOptions.add_argument('-t', "--tree", type=str,
                             help='input tree name', required=True)
requiredOptions.add_argument('-o', '--output', type=str,
                   help="Output directory", required=True)
##############
Options = parser.add_argument_group('Options')
Options.add_argument('-c', '--cont_trait_table', type=str, help="trait table name")
Options.add_argument('-aa', '--aa_align', type=str, help="AA alignment name")
Options.add_argument('-k', '--key_seq', type=str, help="Name of key sequence on which to index the output columns")
Options.add_argument('-tt', '--test_trees', action="store_true", help="Draw test trees to evaluate the discretization scheme")
Options.add_argument('-i', '--low_val_convergent', action="store_true", help="Consider low trait value to be convergent state, i.e. invert the trait")
Options.add_argument('-d', '--det', action="store_true", help="Set to actually run pcoc_det.py")
Options.add_argument('-f', '--float', action="store_true", help="Store trait cutoffs as floating-point values")
Options.add_argument('-p', '--precision', type=float, default=0.0, help="Minimum difference in trait cutoffs between consecutive scenarios to be tested")
Options.add_argument('-hm', '--heatmap', type=str, help="Render heatmap from the latest available set of data and save it here")
#Options.add_argument('-m', '--master_table', action="store_true", help="Print concatenated (master) data table to stdout")
Options.add_argument('-m', '--master_table', type=str, help="Save concatenated (master) data table at...") # TODO: have this go to stdout, but for some reason my stderr goes there too!
Options.add_argument('-mp', '--manhattan', type=str, help="Save Manhattan plot at...")
Options.add_argument('-s', '--sim', action="store_true", help="Set to simulate alignments for the passed tree and continuous trait values")
Options.add_argument('-n', '--aa_noise', action="store_true", help="Add noise to the ancestral profile in simulation")
Options.add_argument('-sd', '--sim_det', action="store_true", help="Set to run shootout detection routine on the simulated alignments")
Options.add_argument('-as', '--analyze_sims', action="store_true", help="Set to calculate true positive and negative rates from the latest set of simulations")
##############

### Option parsing
#global args
args = parser.parse_args()

### Number tree nodes
def init_tree(nf):
    t = Tree(nf)

    #Alternatively we could read a tree from a file into a string "line", and then use:
    # t =  Tree( line )

    nodeId = 0
    for n in t.traverse("postorder"):
        n.add_features(ND=nodeId)
        nodeId = nodeId + 1

    return t


### Perform continuous-trait ancestral state reconstruction. Return list of trait values indexed on node #.
def ancR(tree, tipTraits):

    # Check that the trait list has values for all leaves
    dataError = False
    for leaf in tree.get_leaves():
        if leaf.name not in tipTraits.keys():
            print >> sys.stderr, "No trait data for {}!".format(leaf.name)
            dataError = True
    if dataError:
        print >> sys.stderr, "exiting"
        sys.exit(1)

    # Init list with an element for each node in tree
    nodeTraits = [None] * len(list(tree.traverse()))
    nodeId = 0
    for node in tree.traverse("postorder"): # Amend this to use enumerate()

        if node.name in tipTraits.keys(): # if the node is a tip or its state is specified
            fixTrait = float(tipTraits[node.name])
            node.add_features(trait = fixTrait)
            nodeTraits[nodeId] = fixTrait

        else: # if the state needs to be reconstructed
            #TODO: break out interpolation function. I guess this will require passing the tree by reference?

            weightTraits = dict()
            totalDistance = tree.get_distance(*node.get_descendants()[:2])
            for daughter in node.get_descendants()[:2]:
                weightTraits[daughter.ND] = daughter.trait*(1-(daughter.get_distance(node)/totalDistance))
            recTrait = sum(weightTraits.values())

            node.add_features(trait = recTrait)
            nodeTraits[nodeId] = recTrait

        nodeId = nodeId + 1

    return nodeTraits


def discretize(nodeTraits, precision):
    ### Iteratively discretize continuous trait into binary matrix

    traitMatrix = pd.DataFrame()

    nodeTraitsUnique = np.unique(nodeTraits) # returns sorted unique vals

    # split the trait range into bins that are at least `precision` units apart
    cutoffs = list()
    # this is a result of the tree-fencepost-problem: can have a cutoff between any two values, but then one is the root
    numPossibleCutoffs = len(nodeTraitsUnique) - 1 # -1 or -2?
    for i in range(0, numPossibleCutoffs):
        if nodeTraitsUnique[i+1] - nodeTraitsUnique[i] >= precision: # only if the trait values are far enough apart
            cutoffs.append((nodeTraitsUnique[i]+nodeTraitsUnique[i+1])/2) # place a cutoff in between those values

    # ^REPLACES:
    '''# place cutoffs smack in-between traits in rank order. This does *not* necessarily try a cut on every branch.
    cutoffs = [None] * (len(nodeTraitsUnique)-1)
    for i in range(0, len(cutoffs)):
        cutoffs[i] = int((nodeTraitsUnique[i]+nodeTraitsUnique[i+1])/2)'''

    print >> sys.stderr, "Used precision of " + str(precision) + " trait units to choose " + str(len(cutoffs)) + "/" + str(numPossibleCutoffs) + " possible cutoffs"

    traitMatrix["trait_cont"] = nodeTraits

    for cutoff in np.unique(cutoffs):
        if not args.float: # round the cutoff to int
            opColumn = "trait_cutoff_" + str(int(cutoff))
        else: # keep it floating
            opColumn = "trait_cutoff_" + str(cutoff)
        traitMatrix[opColumn] = pd.cut(traitMatrix["trait_cont"], [min(nodeTraits)]+[cutoff]+[max(nodeTraits)], include_lowest=True, labels=False)

    return traitMatrix

# Consolidate scenarios that are identical under their average header value
def uniqueFilter(traitMatrix):

    #print >> sys.stderr, "pre-dedup-filter:"
    #traitMatrix.to_csv(sys.stderr, sep='\t') #TEST

    keyIndex = 0
    while keyIndex < len(traitMatrix.columns): # need a while loop because this routine is recursive
    #for keyIndex, keyColname in enumerate(list(traitMatrix)):
        keyColname = list(traitMatrix)[keyIndex]
        if len(findall( r'\d+\.*\d*', keyColname)) != 0: # if there is a numeric trait value in the colname
            averageMe = [float(findall( r'\d+\.*\d*', keyColname)[0])] # init a list of names of columns to be consolidated
            keyCol = traitMatrix[keyColname] # this is the column that others will be compared to
            for colname in list(traitMatrix)[keyIndex+1:]: # iterate through the remaining columns
                col = traitMatrix[colname]
                if col is keyCol: # if one matches the key column
                    averageMe.append(float(findall( r'\d+\.*\d*', colname)[0])) # add its name to the list to be averaged
                    traitMatrix = traitMatrix.drop(colname) # and drop it from the matrix.
            # once all duplicates are dropped, rename the key column with the mean cutoff of all the duplicates
            traitMatrix = traitMatrix.rename(index=str, columns={str(keyColname) : "trait_cutoff_"+str(np.mean(averageMe))})
        keyIndex += 1

    #print >> sys.stderr, "post-dedup-filter:"
    #traitMatrix.to_csv(sys.stderr, sep='\t')  #TEST

    return traitMatrix

# Turn the binary trait matrix into a dict of scenario strings for pcoc_det.py
def scenarioStrings(binDf, tree):

    # scenario strings keyed on cutoff value
    scenarios = dict()

    # maybe this could be thinned out by actually making the matrix bool()
    if args.low_val_convergent:
        # invert the trait if specified
        binDf = ~binDf.astype(bool)
    else:
        # otherwise high trait value is considered the convergent state
        binDf = binDf.astype(bool)

    # count and remove scenarios where the root is convergent
    numRootConv = 0
    for col in binDf.columns:
        rootConv = False
        scenario = str()
        for node, state in enumerate(binDf[col]):
            if state:
                if node != binDf.shape[0]-1: #if not the root node
                    parentND = tree.search_nodes(ND=node)[0].up.ND
                    if binDf[col][parentND]:
                        scenario = ','+str(node)+scenario
                    else:
                        scenario = '/' + str(node) + scenario
                else: #if the root node is convergent
                    scenario = ',' + str(node) + scenario
                    rootConv = True
                    numRootConv += 1

        if rootConv == False:
            # key on a sanitized form of the cutoff value
            scenarios[float(findall( r'\d+\.*\d*', col)[0])] = scenario[1:]

    # notify user of outcome of root filter
    print str(numRootConv) + "/" + str(len(scenarios) + numRootConv) + " scenarios eliminated due to convergent root"
    if numRootConv >= (len(scenarios) + numRootConv)/2:
        print "Consider inverting your trait!"

    # to get unique scenarios, I'd have to filter the entries with unique values. But then which key would I use?
    # better to filter unique downstream of the precision filter.
    return scenarios


# take a directory of simulated alignment files produced by pcoc_sim.py and return a list of these filenames
# as Ancestral, Convergent pairs with full paths
def getSimAlignPairs(alignsDir):
    # sorting the filenames is important!
    filenames = sorted(os.listdir(alignsDir))
    # split into duples, which then just need to be ordered
    pairs = [tuple(filenames[i * 2:(i + 1) * 2]) for i in range((len(filenames) + 2 - 1) // 2 )]

    # check each pair to see if it's in the right order
    for i, pair in enumerate(pairs):
        # get the profile numbers for the first member
        aProfile = int(search('A(\d+)', pair[0]).group(1))
        cProfile = int(search('C(\d+)', pair[0]).group(1))
        # if the first member is the convergent alignment
        if aProfile != cProfile:
            # swap the order!
            pairs[i] = pair[::-1]
        # and then add the absolute path
        pairs[i] = tuple([alignsDir + '/' + filename for filename in pairs[i]])

    return pairs


# run a series of pcoc_sim routines on the passed scenario strings and specified tree
def pcocSimulate(scenarios, treedir, simDir, nsites=100, profilePairs=2):

    simAligns = dict()
    s = 1
    for cutoff in sorted(scenarios.keys()):
        try:
            sim_log = open(args.output + "/sim.log", 'a')
        except:
            pass
        # log header
        sim_log.write("~~~ evolutionary scenario " + str(s) + "/" + str(len(scenarios)) + " ~~~\n")
        sim_log.write("trait cutoff = " + str(cutoff) + "\n")
        # command line tracker
        print >> sys.stderr, str(s) + "/" + str(len(scenarios)) + ": " + str(cutoff)
        simArgs = ["pcoc_sim.py",
                   "-td", treedir,
                   "-o", simDir + '/' + str(cutoff).replace('.', '_'),
                   # "-cpu " + str(cpus),
                   "-m", str(scenarios[cutoff]),  # run in manual mode with the generated scenario
                   "-n_sites", str(nsites),
                   "-nb_sampled_couple", str(profilePairs),
                   "--no_clean_seqs"]  # do NOT clean up the simulated alignments!
        # add AA noise if requested
        if args.aa_noise: simArgs.append("--ali_noise")
        subprocess.call(simArgs, stdout=sim_log,
                        stderr=sim_log)  # This will put the scenario string in the main logfile
        sim_log.close()  # give the stream back to parent process (this one)
        s += 1

        #subDir = simDir + '/' + str(cutoff).replace('.', '_')
        #latestRun = sorted(os.listdir(subDir))[-1]  # get the output dir from the latest run
        #alignsDir = subDir + "/" + latestRun + "/Tree_1/sequences/Scenario_1"
        ## get a list of pairs of ancestral/convergent simulated alignments
        #acPairs = getSimAlignPairs(alignsDir)
        #simAligns[cutoff] = acPairs # append the pairs of alignments to dict for return

    # return dict of lists of pairs of alignment filenames, keyed on cutoff value
    #return simAligns


def cpglsDetect(treefile, alignfile, responsefile, keySeq=None, log=sys.stderr):
    # args are parsed within the main() function
    cpglsArgs = ["-t", treefile,
                 "-p", alignfile,
                 "-f", "fasta",
                 "-r", responsefile]
    # "--cpu", "1",
    # "-m"] # calling with plotting option tends to crash
    if keySeq: cpglsArgs = cpglsArgs + ["-k", args.key_seq]

    # return DataFrame()
    return cpgls.main(cpglsArgs, log)


# run a series of pcoc_det routines on the passed scenario strings and specified tree, alignment
def pcocDetect(scenarios, scenDir, treefile=args.tree, alignfile=args.aa_align, threshold=0.8):
    # make a scenarios directory
    try:
        os.mkdir(scenDir)
    except:
        pass

    s = 1
    for cutoff in sorted(scenarios.keys()):
        try:
            det_log = open(args.output + "/det.log", 'a')
        except:
            pass
        # log header
        det_log.write("~~~ evolutionary scenario " + str(s) + "/" + str(len(scenarios)) + " ~~~\n")
        det_log.write("trait cutoff = " + str(cutoff) + "\n")
        # command line tracker
        print >> sys.stderr, str(s) + "/" + str(len(scenarios)) + ": " + str(cutoff)
        detArgs = ["pcoc_det.py",
                   "-t", treefile,
                   "-aa", alignfile,
                   "-o", scenDir + '/' + str(cutoff).replace('.', '_'),
                   # "-cpu " + str(cpus),
                   "-m", str(scenarios[cutoff]),
                   #"--plot_complete_ali",
                   #"--plot",
                   "-f", str(threshold)]
        subprocess.call(detArgs, stdout=det_log,
                        stderr=det_log)  # This will put the scenario string in the main logfile
        det_log.close()  # give the stream back to parent process (this one)
        s = s + 1

    det_log.close()


# gather up all the output files (unfiltered) and put them into a master dataframe
def consolidatePCOCOutput(scenarios, scenDir):
    metaDf = pd.DataFrame()
    propTable = os.path.splitext(os.path.basename(args.aa_align))[:-1][
                    0] + ".results.tsv"  # get the output file from the latest run

    for cutoff in sorted(scenarios.keys()):

        subDir = scenDir + '/' + str(cutoff).replace('.','_')
        latestRun = sorted(os.listdir(subDir))[-1]
        ppFile = subDir + "/" + latestRun + "/" + propTable
        try:
            cutoffPPs = pd.read_table(ppFile)
            cutoffPPs["cutoff"] = cutoff
            metaDf = metaDf.append(cutoffPPs)
        except:
            pass

    return metaDf


def traitTree(traits, mapper):
    ### Take dict of traits and [R,G,B]-returning function
    ### Draw a tree with the continuous trait painted on via a colormapping function

    def rgb2hex(r, g, b):
        hex = "#{:02x}{:02x}{:02x}".format(r, g, b)
        return hex

    for n in tree.traverse():
        if n.is_leaf():
            n.set_style(nstyle_L)
            n.add_face(TextFace(str(n.name)), column=0, position="aligned")
        else:
            n.set_style(nstyle)
        #nd = TextFace(str(n.ND)) # label with node ID
        nd = TextFace(str(int(traits[n.ND]))) # label with rounded continuous trait value

        nd.background.color = rgb2hex(*[int(val) for val in mapper(traits[n.ND], gamma=0.8, scaleMax=255)]) # setup for wl2RGB
        nd.margin_right = 2
        nd.margin_top = 1
        nd.margin_left = 2
        nd.margin_bottom = 1
        nd.border.width = 1
        n.add_face(nd, column=0, position="float")
        n.add_face(TextFace("       "), column=0, position="branch-bottom")

    outfile = args.output + "/test_trees/cont_trait.pdf"
    tree.render(outfile, tree_style=tree_style)
    print >> sys.stderr, outfile
    # no return


def testTrees(scenarios, tree_style):
    ### Draw test trees. This is a modified version of the test routine in pcoc_num_tree.py, stuffed in a for loop
    for cutoff in sorted(scenarios.keys()):
        # this keeps attributes from stacking up in the same tree
        tree = init_tree(args.tree)
        manual_mode_nodes = {}
        manual_mode_nodes = {"T": [], "C": []}
        p_events = scenarios[cutoff].strip().split("/")
        for e in p_events:
            l_e = map(int, e.split(","))
            manual_mode_nodes["T"].append(l_e[0])
            manual_mode_nodes["C"].extend(l_e[1:])

        for n in tree.traverse():
            if n.is_leaf():
                n.set_style(nstyle_L)
                n.add_face(TextFace(str(n.name)), column=0, position="aligned")
            else:
                n.set_style(nstyle)
            nd = TextFace(str(n.ND))

            if manual_mode_nodes:
                if n.ND in manual_mode_nodes["T"]:
                    nd.background.color = "red"
                elif n.ND in manual_mode_nodes["C"]:
                    nd.background.color = "orange"
                else:
                    nd.background.color = "white"
            else:
                nd.background.color = "white"
                nd.background.color = "white"
            nd.margin_right = 2
            nd.margin_top = 1
            nd.margin_left = 2
            nd.margin_bottom = 1
            nd.border.width = 1
            n.add_face(nd, column=0, position="float")
            n.add_face(TextFace("       "), column=0, position="branch-bottom")

        outfile = args.output + "/test_trees/" + str(cutoff).replace('.','_') + ".pdf"
        tree.render(outfile, tree_style=tree_style)
        print >> sys.stderr, outfile
        # no return


def manhattanPlot(df, outPath, keyID=None):
    # at present this is meant to be run once
    import matplotlib.pyplot as plt
    import math

    fig, ax = plt.subplots()
    x = [i+1 for i in range(len(df))]
    y = np.array(df)
    #y = np.array([-math.log10(p) for p in y]) # log-transform
    xlimits = (min(x)-0.5, max(x)+0.5)
    plt.xlim(xlimits)
    plt.xticks(x, x, rotation=90, fontsize=8)

    if keyID:
            ax.set_xlabel("position in " + keyID)
    else:
            ax.set_xlabel("alignment column")
    ax.set_ylabel("total PCOC PP")

    #pThresholds = (1, 0.05, 0.01, 0.001)
    #pThresholds = [-math.log10(p) for p in pThresholds] # log-transform
    pThresholds = (0, 0.8, 0.9, 0.95, 1)
    colors = ["black", "blue", "red", "violet"]
    plt.hlines(pThresholds[:-1], xlimits[0], xlimits[1], colors=colors)

    # for lin
    #masks = [[pThresholds[i] >= el > pThresholds[i + 1] for el in y] for i in range(len(pThresholds) - 1)]
    # for log
    masks = np.array([[pThresholds[i] <= el < pThresholds[i + 1] for el in y] for i in range(len(pThresholds) - 1)])

    for i, mask in enumerate(masks):
            plt.bar(x*mask, y*mask, color=colors[i], width=1)
    
    fig.set_size_inches(30, 5)
    #plt.show()
    plt.savefig(outPath)

    # no return


def key_alignment_columns(alignment, key_seqid):
    '''return list of columns of an alignment that are not gaps in the key sequence'''
    al_length = alignment.get_alignment_length()

    # Get the key sequence
    for seqrec in alignment:
        if seqrec.id == key_seqid:
            keyseq = seqrec.seq

    mappedCols = [None] * len(str(keyseq).replace('-', ''))
    keyindex = 0
    for i in range(al_length):  # only iterate over columns that are not gaps in target seq
        if keyseq[i] != "-":  # meaning anything except gaps
            mappedCols[keyindex] = i + 1  # map the alignment column to the key column. Gotta shift the index!
            keyindex += 1

    # print >> sys.stderr, mappedCols #TEST
    return mappedCols


# Basic tree style
tree_style = TreeStyle()
tree_style.show_leaf_name = False
tree_style.show_branch_length = False
tree_style.draw_guiding_lines = True
tree_style.complete_branch_lines_when_necessary = True

nstyle = NodeStyle()
nstyle["fgcolor"] = "black"
nstyle["size"] = 0

nstyle_L = NodeStyle()
nstyle_L["fgcolor"] = "black"
nstyle_L["size"] = 0

if not os.path.isfile(args.tree):
    print ("%s does not exist" %args.tree)
    sys.exit(1)

tree = init_tree(args.tree)
# not mucking with additive trees yet; ultrametrize the tree and normalize to length 1
tree.convert_to_ultrametric(tree_length=1)

### Load in continuous trait data

def main(wayout):

    # Load in dict of traits keyed on species. Note hardcoding of 'sp' colname!
    tipTraits = pd.read_table(args.cont_trait_table).set_index('sp').transpose().to_dict(orient='r')[0]

    # Generate discrete trait matrix for all cutoff values
    nodeTraits = ancR(tree, tipTraits)
    traitMatrix = discretize(nodeTraits, args.precision)
    traitMatrix = uniqueFilter(traitMatrix) # drop columns as needed to eliminate redundant scenarios
    scenarios = scenarioStrings(traitMatrix.drop(traitMatrix.columns[0], axis=1), tree) # drop first column
    #scenarios = scenarioStrings(traitMatrix.set_index(traitMatrix.columns[0]), tree)  # set first column to index

    #print scenarios #TEST

    simDir = args.output + "/simulations"
    simResultsDir = simDir + "/sim-results"
    scenDir = args.output + "/scenarios"
    if args.test_trees or args.sim or args.det:

        # make a dir for the test trees
        testDir = args.output + "/test_trees"
        try:
            os.mkdir(testDir)
        except:
            pass

        traitTree(nodeTraits, wavelength_to_rgb) # draw trait-colored tree
        testTrees(scenarios, tree_style) # draw test trees

    # if invoking any external PCOC programs, write an ultrametric treefile in a special location
    if args.sim or args.det or args.analyze_sims:
        # pcoc_sim uses a dir full of trees, so I will copy the tree over to its own directory for now.
        # This will be extensible to transformations on the tree (e.g. BL noise) later.
        tdTemp = args.output + "/treetemp"
        try:
            os.mkdir(tdTemp)
        except:
            pass
        # put the (ultrametrized) tree in there
        ultramTree = tdTemp + "/temp.tree"
        tree.write(outfile=ultramTree)

    if args.sim: # if in simulate-alignments mode

        sim_log = open(args.output + "/sim.log", 'a')
        sim_log.write("system args:\n")
        sim_log.write(str(args)) # for the record!
        sim_log.write("continuous trait values:\n")
        pprint(tipTraits, stream=sim_log)  # Log the continuous trait values used

        # make a dir for all that pcoc_sim output
        try:
            os.mkdir(simDir)
        except:
            pass

        print >> sys.stderr, "running C60 sequence simulation"
        # generate simulated alignments and get their filepaths
        pcocSimulate(scenarios, treedir=tdTemp, simDir=simDir, profilePairs=2)

    if args.sim_det:

        simAlignPairs = dict()
        for cutoff in sorted(scenarios.keys()):
            # recover the latest set of simulated sequences
            subDir = simDir + '/' + str(cutoff).replace('.', '_')
            latestRun = sorted(os.listdir(subDir))[-1]  # get the output dir from the latest run
            alignsDir = subDir + "/" + latestRun + "/Tree_1/sequences/Scenario_1"
            # get a list of pairs of ancestral/convergent simulated alignments, in that order
            acPairs = getSimAlignPairs(alignsDir)
            simAlignPairs[cutoff] = acPairs  # append the pairs of alignments to dict for return

        # make dirs for the simulation results output
        try: os.mkdir(simResultsDir)
        except: pass
        try: os.mkdir(simResultsDir + "/PGLS-anc")
        except: pass
        try: os.mkdir(simResultsDir + "/PGLS-con")
        except: pass

        # analyze the simulated alignments: just like below, but on a loop with the simulated alignments
        det_log = open(args.output + "/sim-det.log", 'a')
        det_log.write("system args:\n")
        det_log.write(str(args))
        det_log.write("continuous trait values:\n")
        pprint(tipTraits, stream=det_log)  # Log the continuous trait values used

        for cutoff in sorted(simAlignPairs.keys()):
            for i, pair in enumerate(simAlignPairs[cutoff]):
                # calc CPGLS p-values
                aResults = cpglsDetect(treefile=ultramTree, alignfile=pair[0], responsefile=args.cont_trait_table,
                            keySeq=args.key_seq, log=det_log)
                cResults = cpglsDetect(treefile=ultramTree, alignfile=pair[1], responsefile=args.cont_trait_table,
                            keySeq=args.key_seq, log=det_log)
                # sanitized scenario ID
                saniScen = str(cutoff).replace('.', '_')
                # save the simulation analysis results
                aFile = simResultsDir + "/PGLS-anc/" + saniScen + "-anc" + str(i) + ".tsv"
                aResults.to_csv(aFile, sep='\t')
                cFile = simResultsDir + "/PGLS-con/" + saniScen + "-con" + str(i) + ".tsv"
                cResults.to_csv(cFile, sep='\t')

                # calc PCOC PPs across the entire range of scenarios (fully orthogonal)
                # They get saved in a dir called "sim#_#####" that's the simulation scenario
                # and then a subdir called "#_#####" that's the detection scenario
                pcocDetect(scenarios, treefile=ultramTree, alignfile=pair[0], scenDir=simResultsDir + "/PCOC-anc/sim" + saniScen + "-anc" + str(i))
                pcocDetect(scenarios, treefile=ultramTree, alignfile=pair[1], scenDir=simResultsDir + "/PCOC-con/sim" + saniScen + "-con" + str(i))



    if args.det: # if in determine mode

        det_log = open(args.output + "/det.log", 'a')
        det_log.write("system args:\n")
        det_log.write(str(args))
        det_log.write("continuous trait values:\n")
        pprint(tipTraits, stream=det_log)  # Log the continuous trait values used

        # run CPGLS
        print >> sys.stderr, "running CPGLS detection"
        cpglsResultsDf = cpglsDetect(treefile=ultramTree, alignfile=args.aa_align, responsefile=args.cont_trait_table, keySeq=args.key_seq, log=det_log)

        ## make a dir for all that pcoc_det output
        #try:
        #    os.mkdir(scenDir)
        #except:
        #    pass

        # run pcoc_det
        print >> sys.stderr, "running PCOC detection"
        pcocDetect(scenarios, treefile=ultramTree, scenDir=scenDir)

    if args.heatmap or args.master_table or args.manhattan:
        # get master dataframe of results for all cutoffs
        metaDf = consolidatePCOCOutput(scenarios, scenDir)
        if args.key_seq: # map the PP scores to key sequence positions
            alignment = AlignIO.read(args.aa_align, format="fasta")
            # get key columns of the alignment
            keyCols = key_alignment_columns(alignment, args.key_seq)
            # filter the 2-column data table to contain only these sites
            metaDf = metaDf[metaDf["Sites"].isin(keyCols)]
            # reindex the sites
            metaDf["Sites"] = [keyCols.index(i)+1 for i in metaDf["Sites"]] # note index shift!
        heatmapDf = metaDf.pivot(index="cutoff", columns="Sites", values="PCOC")

    if args.heatmap:
        # Make a heatmap figure. Graphics parameters are all set in `pcoc_cont_heatmap.py`
        rainbow = False #Is the trait in question a visible wavelength to be plotted chromatically?
        heatmapper.heatMapDF(heatmapDf, args.heatmap, rainbow=rainbow)

    if args.manhattan:
        # Make a Manhattan plot

        def sumProbs(probArray):
            # sum columns of non-exclusive probabilities
            # the array casting is important, because numpy uses the operators element-wise
            sums = np.array([0.] * probArray.shape[1])
            for row in probArray:
                sums = sums + np.array(row) - (sums * np.array(row))
            return sums

        #manhattanDf = heatmapDf.sum(axis=0)
        manhattanSeries = sumProbs(np.array(heatmapDf)) # convert to array
        #print sumProbs(np.array(heatmapDf)[:,29:30])
        #manhattanDf /= manhattanDf.max() #normalize
        manhattanPlot(manhattanSeries, args.manhattan, args.key_seq)

    if args.det and args.master_table:
        # transpose the PCOC output
        metaDf = metaDf.pivot(columns="cutoff", index="Sites", values="PCOC")
        # join it with the PGLS output
        metaDf = metaDf.join(cpglsResultsDf)

        # Print master data table. Note that it is transpose of the table sent to heatMapDF()
        metaDf.to_csv(args.master_table, sep='\t')
        #metaDf.pivot(columns="cutoff", index="Sites", values="PCOC").to_csv(args.master_table, sep='\t')

    #print(scenarios) #TEST


if __name__ == "__main__":
        main(sys.stdout)