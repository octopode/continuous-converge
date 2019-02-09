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

##########
# inputs #
##########
startDateTime = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

if __name__ == "__main__":
    
    ### Option defining
    parser = argparse.ArgumentParser(prog="pcoc_cont_scenarios.py", description='')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')

    ##############
    requiredOptions = parser.add_argument_group('Required arguments')
    requiredOptions.add_argument('-t', "--tree", type=str, help='input tree name', required=True)
    requiredOptions.add_argument('-o', '--output', type=str, help="Output directory", required=True)
    requiredOptions.add_argument('-c', '--cont_trait_table', type=str, help="trait table name")
    ##############
    Options = parser.add_argument_group('Options')
    Options.add_argument('-aa', '--aa_align', type=str, help="AA alignment name")
    Options.add_argument('-k', '--key_seq', type=str, help="Name of key sequence on which to index the output columns")
    Options.add_argument('-tt', '--test_trees', action="store_true",
                         help="Draw test trees to evaluate the discretization scheme")
    Options.add_argument('-i', '--invert_trait', action="store_true",
                         help="Invert the binary trait, i.e. assert that low trait value is the convergent state")
    Options.add_argument('-d', '--det', action="store_true", help="Set to actually run pcoc_det.py")
    Options.add_argument('-f', '--float', type=int, default=None, help="Store trait cutoffs as floating-point values")
    Options.add_argument('-p', '--precision', type=float, default=0.0,
                         help="Minimum difference in trait cutoffs between consecutive scenarios to be tested")
    Options.add_argument('-hm', '--heatmap', type=str,
                         help="Render heatmap from the latest available set of data and save it here")
    Options.add_argument('-m', '--master_table', type=str, help="Save collated master data table at...")
    #TODO: have master table go to stdout, but for some reason my stderr goes there too!
    Options.add_argument('-mp', '--manhattan', type=str, help="Save Manhattan plot at...")
    ##############

    ### Option parsing
    #global args
    args = parser.parse_args()

    main(sys.stdout)

### Number tree nodes for consistent reference
def init_tree(nf):
    t = Tree(nf)

    nodeId = 0
    for n in t.traverse("postorder"):
        n.add_features(ND=nodeId)
        nodeId = nodeId + 1

    return t

### Perform BM continuous-trait ancestral state reconstruction. Return list of trait values indexed on node #.
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

### Iteratively discretize continuous trait into binary matrix
def discretize(nodeTraits, precision):

    binDF = pd.DataFrame()

    nodeTraitsUnique = np.unique(nodeTraits) # returns sorted unique vals

    # split the trait range into bins that are at least `precision` units apart
    cutoffs = list()
    numPossibleCutoffs = len(nodeTraitsUnique) - 1
    for i in range(0, numPossibleCutoffs):
        # if the trait values are far enough apart
        if nodeTraitsUnique[i+1] - nodeTraitsUnique[i] >= precision:
            # place a cutoff smack in between those values
            cutoffs.append((nodeTraitsUnique[i]+nodeTraitsUnique[i+1])/2)

    #print >> sys.stderr, "Used precision of " + str(precision) + " trait units to choose " + str(len(cutoffs)) + "/" + str(numPossibleCutoffs) + " possible cutoffs"
    print >> sys.stderr, "# MESSAGE: Used precision of {} trait units to choose {}/{} possible cutoffs".format(precision, len(cutoffs), numPossibleCutoffs)

    # put traits into DF col 0
    binDF["trait_cont"] = nodeTraits

    for cutoff in np.unique(cutoffs):
        if not args.float: # round the cutoff to int
            opColumn = "trait_cutoff_" + str(int(cutoff))
        else: # keep it floating
            opColumn = "trait_cutoff_" + str(cutoff)
        binDF[opColumn] = pd.cut(binDF["trait_cont"], [min(nodeTraits)]+[cutoff]+[max(nodeTraits)], include_lowest=True, labels=False)

    # drop the first column, which contains continuous trait values
    binDF = binDF.drop(binDF.columns[0], axis=1)

    # convert binary values to Boolean
    if args.invert_trait:
        # invert the trait if specified
        binDF = ~binDF.astype(bool)
    else:
        # otherwise high trait value is considered the convergent state
        binDF = binDF.astype(bool)

    return binDF

### Consolidate identical scenarios under their average header value
### This filter is topology-agnostic, so works on the binary trait matrix only
def uniqueColumnFilter(binDF):

    chaff = pd.DataFrame() # init DF for dropped columns
    origNumScenarios = binDF.shape[1]
    keyIndex = 0

    while keyIndex < len(binDF.columns): # need a while loop because this routine is recursive
        keyColName = list(binDF)[keyIndex]
        # if there is a numeric trait value in the colName
        if len(findall( r'\d+\.*\d*', keyColName)) != 0:
            # init a list of names of columns to be consolidated
            averageMe = [float(findall( r'\d+\.*\d*', keyColName)[0])]
            keyCol = binDF[keyColName] # this is the column that others will be compared to
            # iterate through the remaining columns
            for colName in list(binDF)[keyIndex+1:]:
                col = binDF[colName]
                # if one matches the key column
                if col is keyCol:
                    # add its cutoff to the list to be averaged
                    averageMe.append(float(findall( r'\d+\.*\d*', colName)[0]))
                    # and put it in the chaff DF
                    chaff[colName] = binDF[colName]
            # once all duplicates are dropped, rename the key column with the mean cutoff of all the duplicates
            binDF = binDF.rename(index=str, columns={str(keyColName) : "trait_cutoff_"+str(np.mean(averageMe))})
        keyIndex += 1

    numRedundant = chaff.shape[1]
    chaffCutoffs = [float(findall(r'\d+\.*\d*', colName)[0]) for colName in chaff.columns]

    # if any scenarios were non-unique
    if numRedundant:
        # list them
        print >> sys.stderr, "# MESSAGE: {}/{} scenarios were redundant:".format(numRedundant, origNumScenarios)
        print >> sys.stderr, "{}".format(chaffCutoffs)

    return binDF, chaff

### Remove and report on scenarios in which the root is called convergent
### root status is deduced from the row names, which are post-order node IDs
def convergentRootFilter(binDF):

    chaff = pd.DataFrame() # init DF for dropped columns
    origNumScenarios = binDF.shape[1]

    for colNum, colName in enumerate(binDF.columns):
        for node, state in enumerate(binDF[colName]):
            # if node is convergent AND also the root
            if state & (node == binDF.shape[0] - 1):
                # put it in the chaff DF
                #print binDF[colName] #TEST
                chaff[colName] = binDF[colName]
                # and drop it from the matrix
                binDF.drop(colName, inplace=True, axis=1) # Not sure whether drop() can take a whole column or needs a colname!!

    # notify user of outcome of root filter
    numRootConv = chaff.shape[1]
    chaffCutoffs = [float(findall( r'\d+\.*\d*', colName)[0]) for colName in chaff.columns]

    print >> sys.stderr, "# MESSAGE: {}/{} scenarios eliminated due to convergent root:".format(numRootConv, origNumScenarios)
    print >> sys.stderr, "{}".format(chaffCutoffs)

    if numRootConv >= origNumScenarios / 2:
        print >> sys.stderr, "# WARNING: Consider inverting your trait!"

    #print binDF #TEST

    return binDF, chaff

### Use tree and binary trait matrix to generate 'convergent scenario' strings readable by PCOC scripts
def scenarioStrings(binDF, tree):

    # scenario strings keyed on cutoff value
    scenarios = dict()

    for col in binDF.columns:
        scenario = str()
        for node, state in enumerate(binDF[col]):
            # if convergent state is true
            if state & (node < binDF.shape[0] - 1):
                parentND = tree.search_nodes(ND=node)[0].up.ND
                if binDF[col][parentND]:
                    scenario = ',' + str(node) + scenario
                else:
                    scenario = '/' + str(node) + scenario

        # remove leading '/' and add scenario to dict, keyed on cutoff
        scenarios[float(findall(r'\d+\.*\d*', col)[0])] = scenario[1:]

    return scenarios

### Gather up *all* the output files in ScenDir and put them into a master dataframe
def consolidatePCOCOutput(scenarios, scenDir):

    metaDf = pd.DataFrame()
    propTable = os.path.splitext(os.path.basename(args.aa_align))[:-1][
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
def traitTree(traits, mapper, outDir):
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

    #outFile = args.output + "/test_trees/cont_trait.pdf"
    outFile = outDir + "/cont_trait.pdf"
    tree.render(outFile, tree_style=tree_style)
    print >> sys.stderr, outFile
    # no return

### Draw a tree with grey branches and black nodes
def traitTreeMinimal(traits, mapper):
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

    outFile = args.output + "/test_trees/cont_trait.pdf"
    tree.render(outFile, tree_style=ts)
    print >> sys.stderr, outFile
    # no return

### Draw test trees in the standard visual style used by pcoc_num_tree.py
def testTrees(scenarios, outDir):

    ### Draw test trees. This is a modified version of the test routine in pcoc_num_tree.py, stuffed in a for loop
    for cutoff in sorted(scenarios.keys()):
        tree = init_tree(args.tree)
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

        # if --float set, limit number of digits in filename
        if args.float:
            outFile = str(cutoff).replace('.','_')[:np.min([args.float, len(str(cutoff))])] + ".pdf"
        else:
            outFile = str(cutoff).replace('.', '_') + ".pdf"

        # prepend path to filename
        outFile = outDir + '/' + outFile
        tree.render(outFile, tree_style=tree_style)
        print >> sys.stderr, outFile
        # no return

### draw color-coded test trees with node convergent status indicated by colored box
def testTreesMinimal(scenarios, traits, mapper):

    def rgb2hex(r, g, b):
        hex = "#{:02x}{:02x}{:02x}".format(r, g, b)
        return hex

    ### Draw test trees. This is a modified version of the test routine in pcoc_num_tree.py, stuffed in a for loop
    for cutoff in sorted(scenarios.keys()):
        tree = init_tree(args.tree)
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
        outFile = args.output + "/test_trees/" + str(cutoff).replace('.','_')[:np.min([5, len(str(cutoff))])] + ".pdf"
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

# Check treefile
if not os.path.isfile(args.tree):
    print >> sys.stderr, "# ERROR: {} does not exist".format(args.tree)
    sys.exit(1)

# Load treefile
tree = init_tree(args.tree)
# not mucking with additive trees yet; ultrametrize the tree and normalize to length 1
tree.convert_to_ultrametric(tree_length=1)


def main(wayout):

    if args.cont_trait_table is not None:

        ### Generate convergent scenarios

        # Load in dict of traits keyed on species. Note hardcoding of 'sp' colName!
        tipTraits = pd.read_table(args.cont_trait_table).set_index('sp').transpose().to_dict(orient='r')[0]
        # Do BM trait reconstruction
        nodeTraits = ancR(tree, tipTraits)
        # Generate discrete trait matrix for all cutoff values
        binDF = discretize(nodeTraits, args.precision)
        # consolidate redundant scenarios
        binDF = uniqueColumnFilter(binDF)[0] # element [1] is the chaff
        # eliminate convergent root scenarios
        binDF = convergentRootFilter(binDF)[0] # element [1] is the chaff
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
            traitTree(nodeTraits, wavelength_to_rgb, testDir)
            #traitTreeMinimal(nodeTraits, wavelength_to_rgb)  # draw trait-colored tree
            # draw boring test trees that indicate convergent scenario with red/orange
            testTrees(scenarios, testDir)
            #testTreesMinimal(scenarios, nodeTraits, wavelength_to_rgb) # draw test trees with color-coded nodes

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
            metaDf = consolidatePCOCOutput(scenarios, scenDir)
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
