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
import datetime
import numpy as np
import pandas as pd
from ete3 import Tree, NodeStyle, TreeStyle, TextFace, CircleFace
from re import findall
from pprint import pprint
from Bio import AlignIO
from wl2rgb import wavelength_to_rgb
import pcoc_cont_heatmap as heatmapper


### Number tree nodes for consistent reference
def init_tree(nf):
    t = Tree(nf)

    for i, n in enumerate(t.traverse("postorder")):
        n.add_features(ND = i)

    return t


### Perform BM continuous-trait ancestral state reconstruction. Return list of trait values indexed on node #.
def ancR(tree, tipTraits):

    # this command does something funky to the branch lengths!
    #tree.convert_to_ultrametric(tree_length=1)

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
        print >> sys.stderr, "# MESSAGE: Use of {} equal-width bins yielded {} unique trait values:".format(fixNumBins, len(activeBins))

    # with specified bin width, centered on midpoint of range
    if fixBinWidth:
        numBins = int((max(nodeTraits) - min(nodeTraits)) / fixBinWidth) + 1
        midpoint = float(max(nodeTraits) + min(nodeTraits)) / 2.0
        start = midpoint - (float(fixBinWidth) * numBins / 2.0)
        binBounds = [start + (fixBinWidth * i) for i in range(numBins + 1)]
        binBounds[-1] *= 1.001  # include the right boundary
        bIndices = [i-1 for i in np.digitize(nodeTraits, binBounds)]
        activeBins = sorted(list(set(bIndices)))
        print >> sys.stderr, "# MESSAGE: Use of bins {} trait units wide yielded {} unique trait values:".format(fixBinWidth, len(activeBins))

    if fixNumBins or fixBinWidth:
    # now, bin the values and average them
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
            print >> sys.stderr, "# {} <- {}".format(avg, init)

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
            for colName, series in binDF.iteritems():
                if tuple(series) == col:
                    # append-in-place no good here
                    toAverage[i] = toAverage[i] + [colName]

        # average each list in the list of lists
        avgCutoffs = [np.mean(cuts) for cuts in toAverage]

        # list them
        print >> sys.stderr, "# MESSAGE: {} scenarios were consolidated into {} unique scenarios:".format(numInitCols,
                                                                                                          numUniqueCols)
        for avg, init in zip(avgCutoffs, toAverage):
            print >> sys.stderr, "# {} <- {}".format(avg, init)

        # construct consolidated dataframe
        binDF = pd.DataFrame(columns = avgCutoffs, data = zip(*uniqueCols))
        # sort the new averaged columns
        binDF.sort_index(axis = 1, inplace = True)

    else:

        print >> sys.stderr, "# MESSAGE: Scenarios for all cutoffs are unique"

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
    print >> sys.stderr, "# MESSAGE: {}/{} scenarios eliminated due to convergent root:".format(len(convergentCutoffs), len(origCutoffs))
    print >> sys.stderr, "{}".format(convergentCutoffs)

    if len(convergentCutoffs) >= len(origCutoffs) / 2:
        print >> sys.stderr, "# WARNING: Consider inverting your trait!"

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

    return scenarios


### Gather up *all* the output files in ScenDir and put them into a master dataframe
def consolidatePCOCOutput(scenarios, scenDir, aa_align):

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
            print >> sys.stdout, "# WARNING: Cutoff " + str(cutoff) + " not loaded"

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

## ETE3 TREE-VIZ FUNCTIONS ##

# basic tree style
tree_style = TreeStyle()
tree_style.show_leaf_name = False
tree_style.show_branch_length = False
tree_style.draw_guiding_lines = True
tree_style.complete_branch_lines_when_necessary = True

# make tree grow upward
tree_style.rotation = 270
# and make it appear ultrametric (which it is!)
tree_style.optimal_scale_level = "full"

# internal node style
nstyle = NodeStyle()
nstyle["fgcolor"] = "black"
nstyle["size"] = 0

# terminal node style
nstyle_L = NodeStyle()
nstyle["fgcolor"] = "black"
nstyle["size"] = 0

### Draw a tree with nodes color-coded and labeled by trait value
def traitTree(tree, traits, mapper, outDir, float=0):
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
        if round == 0:
            nd = TextFace(str(round(traits[n.ND]))) # label with rounded continuous trait value
        else:
            nd = TextFace(str(round(traits[n.ND], float)))  # label with rounded continuous trait value

        nd.background.color = rgb2hex(*[int(val) for val in mapper(traits[n.ND], gamma=0.8, scaleMax=255)]) # setup for wl2RGB
        nd.margin_right = 2
        nd.margin_top = 1
        nd.margin_left = 2
        nd.margin_bottom = 1
        nd.border.width = 1
        n.add_face(nd, column=0, position="float")
        n.add_face(TextFace("       "), column=0, position="branch-bottom")

    outFile = outDir + "/cont_trait.pdf"
    tree.render(outFile, tree_style=tree_style)
    print >> sys.stderr, outFile
    # no return

### Draw a tree with grey branches and black nodes
def traitTreeMinimal(tree, traits, mapper, output):
    ### Take dict of traits and [R,G,B]-returning function
    ### Draw a tree with the continuous trait painted on via a colormapping function
    def rgb2hex(r, g, b):
        hex = "#{:02x}{:02x}{:02x}".format(r, g, b)
        return hex

    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.rotation = 270
    ts.complete_branch_lines_when_necessary = False
    #ts.optimal_scale_level = "full"
    ts.scale = 800

    # default NodeStyle
    nstyle = NodeStyle()
    nstyle["size"] = 0
    nstyle["hz_line_color"] = "grey"
    nstyle["vt_line_color"] = "grey"
    nstyle["hz_line_width"] = 3
    nstyle["vt_line_width"] = 3

    for n in tree.traverse():
        chroma = rgb2hex(*[int(val) for val in mapper(traits[n.ND], gamma=0.8, scaleMax=255)])  # setup for wl2RGB
        nf = CircleFace(radius = 10, color = 'none', style='circle', label=None)
        n.set_style(nstyle)
        n.add_face(nf, column=0, position='branch-top')

    outFile = output + "/test_trees/cont_trait.pdf"
    tree.render(outFile, tree_style=ts)
    print >> sys.stderr, outFile
    # no return

### Draw test trees in the standard visual style used by pcoc_num_tree.py
def testTrees(tree, scenarios, outDir, treePath, floatSwitch=0):

    ### Draw test trees. This is a modified version of the test routine in pcoc_num_tree.py, stuffed in a for loop
    for cutoff in sorted(scenarios.keys()):
        tree = init_tree(treePath)
        # not mucking with additive trees yet; ultrametrize the tree and normalize to length 1
        tree.convert_to_ultrametric(tree_length=1)
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

        outFile = str(round(cutoff, floatSwitch)).replace('.','_') + ".pdf"

        # prepend path to filename
        outFile = outDir + '/' + outFile
        tree.render(outFile, tree_style=tree_style)
        print >> sys.stderr, outFile
        # no return

### draw color-coded test trees with node convergent status indicated by colored box
def testTreesMinimal(tree, scenarios, traits, mapper, treePath, output, floatSwitch = 0):

    def rgb2hex(r, g, b):
        hex = "#{:02x}{:02x}{:02x}".format(r, g, b)
        return hex

    ### Draw test trees. This is a modified version of the test routine in pcoc_num_tree.py, stuffed in a for loop
    for cutoff in sorted(scenarios.keys()):
        tree = init_tree(treePath)
        # not mucking with additive trees yet; ultrametrize the tree and normalize to length 1
        tree.convert_to_ultrametric(tree_length=1)

        # read scenario into a dict
        manual_mode_nodes = {"T": [], "C": []}
        p_events = scenarios[cutoff].strip().split("/")
        for e in p_events:
            l_e = map(int, e.split(","))
            manual_mode_nodes["T"].append(l_e[0])
            manual_mode_nodes["C"].extend(l_e[1:])

        ts = TreeStyle()
        # ts.allow_face_overlap = True
        ts.show_leaf_name = False
        ts.rotation = 270
        ts.complete_branch_lines_when_necessary = False
        # ts.optimal_scale_level = "full"
        ts.scale = 800

        for n in tree.traverse():

            # default NodeStyle
            nstyle = NodeStyle()
            nstyle["size"] = 0
            nstyle["hz_line_color"] = "none"
            nstyle["vt_line_color"] = "none"
            nstyle["hz_line_width"] = 3
            nstyle["vt_line_width"] = 3

            # colored faces
            chroma = rgb2hex(*[int(val) for val in mapper(traits[n.ND], gamma=0.8, scaleMax=255)])  # setup for wl2RGB
            nf = CircleFace(radius=10, color=chroma, style='circle', label=None)

            # scenario-dependent features
            if manual_mode_nodes:
                # if transition node
                if n.ND in manual_mode_nodes["T"]:
                    #nstyle["hz_line_color"] = "orange"
                    nf.inner_border.width = 4
                    nf.inner_border.color = 'red'
                # if convergent node
                elif n.ND in manual_mode_nodes["C"]:
                    #nstyle["hz_line_color"] = "violet"
                    nf.inner_border.width = 4
                    nf.inner_border.color = 'white'
                # if ancestral
                else:
                    nstyle["hz_line_color"] = "none"

                n.set_style(nstyle)
                n.add_face(nf, column=0, position='branch-top')

        # limiting number of digits
        outFile = output + "/test_trees/" + str(round(cutoff, floatSwitch)).replace('.','_') + ".pdf"
        tree.render(outFile, tree_style=ts)
        print >> sys.stderr, outFile
        # no return

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


def main(wayout):

    ##########
    # inputs #
    ##########
    startDateTime = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

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
    Options.add_argument('-d', '--det', action="store_true", help="Set to actually run pcoc_det.py")
    Options.add_argument('-f', '--float', type=int, default=0, help="Round traits to f decimals in filenames and figures")
    Options.add_argument('-nb', '--num_bins', type=int, default=0, help="Average continuous trait into n equal-width bins, 0 = no binning")
    Options.add_argument('-bw', '--bin_width', type=float, default=0, help="Average continuous trait into bins of width w, 0 = no binning")
    #Options.add_argument('-p', '--precision', type=float, default=0.0,
    #                     help="Minimum difference in trait cutoffs between consecutive scenarios to be tested")
    Options.add_argument('-hm', '--heatmap', type=str,
                         help="Render heatmap from the latest available set of data and save it here")
    Options.add_argument('-m', '--master_table', type=str, help="Save collated master data table at...")
    # TODO: have master table go to stdout, but for some reason my stderr goes there too!
    Options.add_argument('-mp', '--manhattan', type=str, help="Save Manhattan plot at...")
    ##############

    ### Option parsing
    # global args
    args = parser.parse_args()
    
    # Check treefile
    if not os.path.isfile(args.tree):
        print >> sys.stderr, "# ERROR: {} does not exist".format(args.tree)
        sys.exit(1)
    
    # Load treefile
    tree = init_tree(args.tree)
    # not mucking with additive trees yet; ultrametrize the tree and normalize to length 1
    tree.convert_to_ultrametric(tree_length=1)
    
    if args.cont_trait_table is not None:

        ### Generate convergent scenarios

        # Load in dict of traits keyed on species. Note hardcoding of 'sp' colName!
        tipTraits = pd.read_table(args.cont_trait_table).set_index('sp').transpose().to_dict(orient='r')[0]
        # Do BM trait reconstruction
        nodeTraits = ancR(tree, tipTraits)
        # get cutoffs using specified binning
        cutoffs = binBy(nodeTraits, fixNumBins = args.num_bins, fixBinWidth = args.bin_width)[0]
        # Generate discrete trait matrix for all cutoff values
        binDF = discretize(nodeTraits, cutoffs, args.invert_trait)
        # consolidate redundant scenarios
        #binDF = uniqueColumnFilter(binDF)[0] # element [1] is the chaff
        binDF = uniqueScenarios(binDF)
        # eliminate convergent root scenarios
        #binDF = convergentRootFilter(binDF)[0] # element [1] is the chaff
        binDF = convergentRootFilter(binDF)
        # convert binary dataframe into scenario strings
        scenarios = scenarioStrings(binDF, tree)

        scenDir = args.output + "/scenarios"

        if args.test_trees or args.det:

            # make a dir for the test trees
            testDir = args.output + "/test_trees"
            try:
                os.mkdir(testDir)
            except:
                pass

            print >> sys.stderr, "# MESSAGE: Drawing trees..."
            # draw trait-colored tree
            traitTree(tree, nodeTraits, wavelength_to_rgb, testDir)
            #traitTreeMinimal(tree, nodeTraits, wavelength_to_rgb, args.output)  # draw trait-colored tree
            # draw boring test trees that indicate convergent scenario with red/orange
            testTrees(tree, scenarios, testDir, args.tree, args.float)
            #testTreesMinimal(tree, scenarios, nodeTraits, wavelength_to_rgb, args.tree, args.output) # draw test trees with color-coded nodes

        # if in determine mode
        if args.det:
            # pcoc_det requires a threshold argument, but this script doesn't use the filtered results
            # so this value doesn't matter
            threshold = 0.8
            # make a dir for all that pcoc_det output
            try:
                os.mkdir(scenDir)
            except:
                pass

            # open the pcoc_det logfile and tell user about it
            print >> sys.stderr, "# MESSAGE: Opening site detection logfile:"
            detLogPath = args.output + "/det_{}.log".format(startDateTime)
            detLogHandle = open(detLogPath, 'a')
            print >> sys.stderr, detLogPath

            # log the continuous trait values used
            detLogHandle.write("continuous trait values:\n")
            pprint(tipTraits, stream=detLogHandle)

            # run pcoc_det.py for each convergent scenario
            print >> sys.stderr, "# MESSAGE: Detecting convergent sites for scenarios..."
            for num, cutoff in enumerate(sorted(scenarios.keys())):
                try: detLogHandle = open(detLogPath, 'a')
                except: pass
                detLogHandle.write("~~~ convergent scenario {}/{} ~~~\n".format(num+1, len(scenarios)))
                detLogHandle.write("trait cutoff = {}\n".format(cutoff))
                # if I don't reopen the file here, the above entries end up below the pcoc_det run. #TODO follow this up
                detLogHandle.close()
                detLogHandle = open(detLogPath, 'a')
                detArgs = ["pcoc_det.py",
                           "-t", args.tree,
                           "-aa", args.aa_align,
                           "-o", scenDir + '/' + str(cutoff).replace('.','_'),
                           "-m", str(scenarios[cutoff]),
                           #"--plot_complete_ali",
                           #"--plot",
                           "-f", str(threshold)]
                print >> sys.stderr, "{}/{}: {}".format(num+1, len(scenarios), cutoff)
                subprocess.call(detArgs, stdout=detLogHandle, stderr=detLogHandle)
                # give the stream back to parent process (this one)
                detLogHandle.close()

            detLogHandle.close()

        ### Collate, visualize results

        if args.heatmap or args.master_table or args.manhattan:
            # get master dataframe of results for all cutoffs
            metaDf = consolidatePCOCOutput(scenarios, scenDir, args.aa_align)
            if args.key_seq: # map the PP scores to key sequence positions
                alignment = AlignIO.read(args.aa_align, format="fasta")
                # get key columns of the alignment
                keyCols = keyAlignmentColumns(alignment, args.key_seq)
                # filter the 2-column data table to contain only these sites
                metaDf = metaDf[metaDf["Sites"].isin(keyCols)]
                # reindex the sites
                metaDf["Sites"] = [keyCols.index(i)+1 for i in metaDf["Sites"]] # note index shift!
                heatmapDf = metaDf.pivot(index="cutoff", columns="Sites", values="PCOC")

        if args.heatmap:
            print >> sys.stderr, "# MESSAGE: Drawing heatmap:"
            # Make a heatmap figure. Graphics parameters are all set in `pcoc_cont_heatmap.py`
            heatmapPath = args.output + '/' + args.heatmap
            rainbow = True #Is the trait in question a visible wavelength to be plotted chromatically?
            heatmapper.heatMapDF(heatmapDf, heatmapPath, rainbow=rainbow)
            # tell user where it was put
            print >> sys.stderr, heatmapPath

        if args.manhattan:

            ### Sum columns of non-exclusive probabilities
            def sumProbs(probArray):
                # the array casting is important, because numpy uses the operators element-wise
                sums = np.array([0.] * probArray.shape[1])
                for row in probArray:
                    sums = sums + np.array(row) - (sums * np.array(row))
                return sums

            print >> sys.stderr, "# MESSAGE: Drawing Manhattan plot:"
            manhattanSeries = sumProbs(np.array(heatmapDf)) # convert to array
            manhattanPath = args.output + '/' + args.manhattan
            manhattanPlot(manhattanSeries, manhattanPath, args.key_seq)
            # tell user where it was put
            print >> sys.stderr, manhattanPath

        if args.master_table:
            print >> sys.stderr, "# MESSAGE: Saving master table of results:"
            # Print master data table. Note that it is transpose of the table sent to heatMapDF()
            metaDf.pivot(columns="cutoff", index="Sites", values="PCOC").to_csv(args.master_table, sep='\t')
            # tell user where it was put
            print >> sys.stderr, args.master_table

if __name__ == "__main__":
    main(sys.stdout)