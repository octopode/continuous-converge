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
from ete3 import Tree, NodeStyle, TreeStyle, TextFace, CircleFace
from matplotlib.colors import LinearSegmentedColormap
from re import findall
from pprint import pprint
from Bio import AlignIO
from wl2rgb import wavelength_to_rgb
import pcoc_cont_heatmap as heatmapper
from time import sleep

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
    numPossibleCutoffs = len(nodeTraitsUnique) - 1
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


def uniqueFilter(traitMatrix):
    # Consolidate scenarios that are identical under their average header value

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


def consolidatePCOCOutput(scenarios, scenDir):
    # gather up all the output files (unfiltered) and put them into a master dataframe
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
            print >> sys.stdout, "Warning: cutoff " + str(cutoff) + " not loaded"

    return metaDf

def traitTreeOld(traits, mapper):
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


def traitTree(traits, mapper):
    ### Take dict of traits and [R,G,B]-returning function
    ### Draw a tree with the continuous trait painted on via a colormapping function
    def rgb2hex(r, g, b):
        hex = "#{:02x}{:02x}{:02x}".format(r, g, b)
        return hex

    ts = TreeStyle()
    #ts.allow_face_overlap = True
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

    #for n in tree.traverse():
    #    chroma = rgb2hex(*[int(val) for val in mapper(traits[n.ND], gamma=0.8, scaleMax=255)]) # setup for wl2RGB
    #    nd = CircleFace(radius=3, color=chroma)  # label with rounded continuous trait value
    #    nd.hz_align = 1
    #    #add_face_to_node(nd, n, column=0, aligned=False, position='branch-right')
    #    n.add_face(nd, column=0, position='float')

    outfile = args.output + "/test_trees/cont_trait.pdf"
    tree.render(outfile, tree_style=ts)
    print >> sys.stderr, outfile
    # no return


def testTrees(scenarios, traits, mapper):

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
        outfile = args.output + "/test_trees/" + str(cutoff).replace('.','_')[:np.min([5, len(str(cutoff))])] + ".pdf"
        tree.render(outfile, tree_style=ts)
        print >> sys.stderr, outfile
        # no return



def testTreesOld(scenarios, tree_style):
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

        # limiting number of digits
        outfile = args.output + "/test_trees/" + str(cutoff).replace('.','_')[:np.min([5, len(str(cutoff))])] + ".pdf"
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
    plt.xticks(x, x, rotation=90, fontsize=7)
    plt.yticks(fontsize=20)

    if keyID:
            ax.set_xlabel("position in " + keyID, fontsize=20)
    else:
            ax.set_xlabel("alignment column")
    ax.set_ylabel("total PCOC PP", fontsize=20)

    #pThresholds = (1, 0.05, 0.01, 0.001)
    #pThresholds = [-math.log10(p) for p in pThresholds] # log-transform
    #pThresholds = (0, 0.8, 0.9, 0.95, 1)
    #colors = ["black", "blue", "red", "violet"]
    pThresholds = (0, 0.8, 1)
    colors = ["white", (0,1,0.05)]
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

# make tree grow upward
tree_style.rotation = 270
# and make it appear ultrametric (which it is!)
tree_style.optimal_scale_level = "full"

nstyle = NodeStyle()
nstyle["fgcolor"] = "black"
nstyle["size"] = 0

nstyle_L = NodeStyle()
nstyle["fgcolor"] = "black"
nstyle["size"] = 0

if not os.path.isfile(args.tree):
    print ("%s does not exist" %args.tree)
    sys.exit(1)

tree = init_tree(args.tree)
# not mucking with additive trees yet; ultrametrize the tree and normalize to length 1
tree.convert_to_ultrametric(tree_length=1)

### Load in continuous trait data

def main(wayout):
    if args.cont_trait_table is not None:
        # Load in dict of traits keyed on species. Note hardcoding of 'sp' colname!
        tipTraits = pd.read_table(args.cont_trait_table).set_index('sp').transpose().to_dict(orient='r')[0]

        # Generate discrete trait matrix for all cutoff values
        nodeTraits = ancR(tree, tipTraits)
        traitMatrix = discretize(nodeTraits, args.precision)
        traitMatrix = uniqueFilter(traitMatrix) # drop columns as needed to eliminate redundant scenarios
        scenarios = scenarioStrings(traitMatrix.drop(traitMatrix.columns[0], axis=1), tree) # drop first column

        #print scenarios #TEST

        scenDir = args.output + "/scenarios"
        if args.test_trees or args.det:

            # make a dir for the test trees
            testDir = args.output + "/test_trees"
            try:
                os.mkdir(testDir)
            except:
                pass

            traitTree(nodeTraits, wavelength_to_rgb) # draw trait-colored tree
            testTrees(scenarios, nodeTraits, wavelength_to_rgb) # draw test trees

        if args.det: # if in determine mode
            threshold = 0.8 #test
            # make a dir for all that pcoc_det output
            try:
                os.mkdir(scenDir)
            except:
                pass

            det_log = open(args.output + "/det.log", 'a')
            det_log.write("continuous trait values:\n")
            pprint(tipTraits, stream=det_log)  # Log the continuous trait values used
            s = 1
            for cutoff in sorted(scenarios.keys()):
                try: det_log = open(args.output + "/det.log", 'a')
                except: pass
                detArgs = ["pcoc_det.py",
                           "-t", args.tree,
                           "-aa", args.aa_align,
                           "-o", scenDir + '/' + str(cutoff).replace('.','_'),
                           # "-cpu " + str(cpus),
                           "-m", str(scenarios[cutoff]),
                           "--plot_complete_ali",
                           "--plot",
                           "-f", str(threshold)]
                det_log.write("~~~ evolutionary scenario " + str(s) + "/" + str(len(scenarios)) + " ~~~\n")
                det_log.write("trait cutoff = " + str(cutoff) + "\n")
                sleep(0.05)
                print str(s) + "/" + str(len(scenarios)) + ": " + str(cutoff)
                subprocess.call(detArgs, stdout=det_log, stderr=det_log) # This will put the scenario string in the main logfile
                det_log.close() # give the stream back to parent process (this one)
                s = s+1

            det_log.close()

        # all below is for data collation, visualization

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

        if args.master_table:
            # Print master data table. Note that it is transpose of the table sent to heatMapDF()
            metaDf.pivot(columns="cutoff", index="Sites", values="PCOC").to_csv(args.master_table, sep='\t')

        #print(scenarios) #TEST


if __name__ == "__main__":
        main(sys.stdout)