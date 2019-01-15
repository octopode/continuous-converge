###!/usr/bin/env python

"""pgls.py
v2.3 (c)2018 Jacob Winnikoff
jwinnikoff@mbari.org

Implements Phylogenetically Generalized Least Squares (PGLS) regression as a standalone utility or thru portable methods.

Input: Newick tree and either 2 response tables or one response table and one FASTA alignment
Output: tab-delimited table of p-values for different response correlations

Usage:

To run CPGLS (column-PGLS on a protein alignment, an experimental method designed by the author):
pgls2_2.py -r trait/traittable.tsv -p seq/prot_algt.fasta -t tree/gene_or_species.tree > prot_trait_pValues.tsv

^Note that CPGLS requires a substitution rate matrix appropriate to your MSA; pass with -s, default is BLOSUM62.

To run regular PGLS using a table of quantitative traits as predictors:
pgls2_2.py -r trait/dependent_traittable.tsv -p predictor_traittable.tsv -f table -t tree/gene_or_species.tree > predictor_dependent_pValues.tsv

^-f table tells the program to expect a predictor in TSV format and switches the program into quantitative-predictor mode.

This program can also draw Manhattan plots and plot phylogenetically corrected statistical distributions/tests.
See -b, -m arguments under pgls.py -h.
"""

import sys
import os
import argparse
import csv
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
import statsmodels.stats.api as sms
import scipy.stats as st
import matplotlib.pyplot as plt
import warnings
from ete3 import Tree#, NodeStyle, TreeStyle, TextFace
from Bio import AlignIO
from multiprocessing import pool, cpu_count
from math import factorial
from tqdm import tqdm
from time import sleep

# other files in project
import blosum # import the whole thing so the matrix dict will load

## global declarations
aalphabet = list('ARNDCQEGHILKMFPSTWYV')  # init our AA alphabet. Order must match with substitution matrix!

### take an IUPAC AA code, return 2O-D binary vector
def aa2vec(aa, aalphabet):
	#return [0b1 if x in aa else -0b1 for x in aalphabet] # [-1, 1] coding
	return [0b1 if x in aa else 0b0 for x in aalphabet] # [0, 1] coding


### take a numpy array and make every row sum to 1
def normalizeRowWise(matrix, norm=1):
	return np.array([row * (norm / sum(row)) for row in matrix])


### what it looks like!
def nCr(n,r):
	f = factorial
	return f(n) / f(r) / f(n - r)


### Take 2 sets of Cartesian points along the same domain and return the x-coordinates of their intersections
### note that presently, this rounds x-values in the negative direction
def solveIntersect(xs, ys1, ys2):

	if min(len(ys1), len(ys2)) < len(xs):
		print >> sys.stderr, "ERR: y-values do not span domain!"
		return 1

	ys1 = np.array(ys1)
	ys2 = np.array(ys2)
	ydiff = ys2 - ys1

	# get the indices that proceed sign changes in the array of differences
	crossesZeroAt = np.where(np.diff(np.sign(ydiff)))[0]
	# get the corresponding x-vals
	intersectXs = [xs[i] for i in crossesZeroAt]
	# and return them
	return intersectXs


### Take a tree object and Pagel's lambda, return a phylogenetic variance-covariance matrix and the axis order
### one can also pass a list of taxa, which allows (1) use of a subset of taxa in the tree, (2) custom order of taxa along the matrix axes
# TODO: implement more up-to-date phylogenetic signal strength parameters
def tree2covMatrix(tree, taxa=None, pagel=1):

	# get root node
	root = tree.get_tree_root()

	# init list of leaves (treeNodes with hex IDs)
	leaves = list()
	if len(taxa): # if a list of taxa is passed
		for taxon in taxa:
			# get the node with that name
			node = tree.search_nodes(name=taxon) # note that this method always returns a list!
			# check to make sure the taxon name is unique in the tree
			if len(node) == 0:
				print >> sys.stderr, "# WARNING: Taxon \'{}\' not found in tree, will be dropped from analysis".format(taxon)
			elif len(node) > 1:
				print >> sys.stderr, "# ERROR: {} nodes in the tree are named \'{}\', exiting".format(len(node), taxon)
				sys.exit(1)
			else: # if search returns exactly 1 node as it should:
				if node[0].is_leaf(): leaves.append(node[0])

	else: # if no list of taxa is passed, use all the leaves on the tree in post-order
		for node in tree.traverse("postorder"): # traverse order is important, or the root will be misidentified!
			if node.is_leaf(): leaves.append(node)

	# init covariance matrix, square with number of leaves as dimensions
	covMatrix = np.empty((len(leaves), len(leaves)))

	for i, inode in enumerate(leaves): # row major
		for j, jnode in enumerate(leaves): # then iterate columns
			covMatrix[i, j] = root.get_distance(tree.get_common_ancestor([inode, jnode])) # enter distance between root and LCA

	# apply Pagel's lambda transformation. It's biological rubbish, but simple.
	covMatrix = covMatrix.dot(pagel * np.linalg.inv(covMatrix))
	# return the lambda-transformed variance-covariance matrix and also the ordered list of taxon names
	return covMatrix, [leaf.name for leaf in leaves]


### Perform continuous-trait ancestral state reconstruction. Return list of trait values indexed on node #.
def ancR(tree, tipTraits):

	# initialize the tree with durable node ID
    nodeId = 0
    for n in tree.traverse("postorder"):
        n.add_features(ND=nodeId)
        nodeId = nodeId + 1

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


### Take variance-covariance matrix, column of AAs, substitution rate matrix and weighting coefficient s
def adjustCovMatrix(C, D, S, s=1):
	# C is the "PHYLOGENETIC covariance matrix"
	# D is a "DAYHOFF-like substitution rate matrix"
	# S is the exogenous "STATE matrix", n rows x 20 cols

	B = S.dot(D).dot(S.T)
	#B = S.dot(D).dot(S.T) + S.dot(D.T).dot(S.T) # this "BIOCHEMICAL transition probability matrix" is n x n, to be added elementwise with C
	B = s * B # apply a weighting parameter to the biochemical covariance matrix.
	# s = 0 -> no difference in AA substitution rates; s = 1 -> maximal difference in AA substitution rates as defined in D

	# return an element-wise non-exclusive sum of the covariance matrices
	return B + C - (B * C) # i.e. p + q - pq

### Load a trait table from a CSV. If desired, return the table with only specified taxa in a specified order
def makeNiceTraitTable(fname, orderedTaxa=None):
	# Load the quantitative predictor values into a DataFrame and sort it to match the order of the cov matrix
	traitTable = pd.read_table(fname, sep='\t')

	if orderedTaxa: # if an ordered list of taxa is provided
		# find entries for taxa that are not passed, and drop them to avoid downstream trouble
		traitTable = traitTable.drop(traitTable[~traitTable[traitTable.columns[0]].isin(orderedTaxa)].index)
		# set column 0 (the taxon ID column) to categorical (aka a factor) and order the levels per orderedLeaves
		traitTable[traitTable.columns[0]] = pd.Categorical(traitTable[traitTable.columns[0]], orderedTaxa)
		# sort the table in the order of orderedLeaves,

	# make taxon the index
	traitTable.set_index(traitTable.columns[0], inplace=True)

	return traitTable

### do a PGLS regression for a column of amino acids, corrected using the passed substitution rate matrix
### return p-value of the regression slope, as well as the response value for the optimal AA state break (calculated by minimal set overlap)
def cpgls(responseList, intercept, phyCovMatrix, subMatrix, colVector, subWeight=1, colnum=0):
	# colnum is dummy for now

	exog = np.array([aa2vec(aa, aalphabet) for aa in colVector.flatten().tolist()])  # exogenous "STATE matrix": n rows x 20 cols
	covMatrix = adjustCovMatrix(phyCovMatrix, subMatrix, exog, subWeight)

	# using the formula API
	data = pd.DataFrame()
	#data["response"] = responseList # forces intercept thru 0
	# fix the intercept at the passed value
	data["response"] = [response - intercept for response in responseList]
	data["aa"] = colVector.flatten()
	return smf.gls(formula=("response ~ aa + 0"), data=data, sigma=covMatrix)

### Do a regular old PGLS for one quantitative trait predicted by another quantitative trait
def pgls(responseList, predictor, phyCovMatrix, intercept=True):
	#predictor = np.array(predictor).T # transpose it to make sure sm doesn't choke
	# TODO: should I set an intercept, or no?
	# Add a column of 1s at the beginning of exog so the model incorporates an intercept
	if intercept:
		predictor = np.array([(1, row) for row in predictor])  # now it is n x 2
	# return the model
	return sm.GLS(endog=responseList, exog=predictor, sigma=phyCovMatrix)

### Return list of columns of an alignment that are not gaps in the key sequence
def key_alignment_columns(alignment, key_seqid):

    algtLength = alignment.get_alignment_length()

    # Get the key sequence
    for seqrec in alignment:
        if seqrec.id == key_seqid:
            keyseq = seqrec.seq

    mappedCols = [None] * len(str(keyseq).replace('-', ''))
    keyindex = 0
    for i in range(algtLength):  # only iterate over columns that are not gaps in target seq
        if keyseq[i] != "-":  # meaning anything except gaps
            mappedCols[keyindex] = i + 1  # map the alignment column to the key column. Gotta shift the index!
            keyindex += 1

    # print >> sys.stderr, mappedCols #TEST
    return mappedCols

### Take a list of indices, a list of scores and a tuple of cutoffs and print a Manhattan plot to outPath
### this plotting function IS threadsafe AFAICT
def manhattanPlotBlackOnWhite(x, y, outPath, pThresholds=(0.05, 0.01, 0.001), keyID=None):

	xNames = x
	x = range(len(x)) # so you can have a non-numeric x-axis, but they still get plotted in order
	fig, ax = plt.subplots()
	y = [-np.log10(p) for p in y] # log-transform p-values
	xlimits = (min(x)-0.5, max(x)+0.5)
	ax.set_xlim(xlimits)
	ax.set_xticks(x, minor=True)
	# hide the default tick labels; this is stupid
	ax.set_xticklabels(x, visible=False)
	# show the correct tick labels
	ax.set_xticklabels(xNames, minor=True, rotation=90, fontsize=8)

	if keyID:
		ax.set_xlabel("position in " + keyID, fontsize=20)
	else:
		ax.set_xlabel("predictor column", fontsize=20)
	ax.set_ylabel("-log10(P)", fontsize=20)

	# add upper and lower thresholds so those bars don't disappear!
	pThresholds = [1] + pThresholds + (1e-100,)
	pThresholds = [-np.log10(p) for p in pThresholds] # log-transform alpha-values

	# set the threshold colors
	#colors = ["black", "blue", "red", "violet"] # boring colors
	#colors = ["black", "5CBF9F", "1B56CC", "A0242A"] # mpl doesn't take hex!
	colors = ["black", (0.3608, 0.7490, 0.6235), (0.1059, 0.3373, 0.8000), (0.6275, 0.1412, 0.1647)]# DEEPC colors!
	# and draw threshold lines
	ax.hlines(pThresholds[:-1], xlimits[0], xlimits[1], colors=colors)

	# for lin
	#masks = [[pThresholds[i] >= el > pThresholds[i + 1] for el in y] for i in range(len(pThresholds) - 1)]
	# for log
	masks = np.array([[pThresholds[i] <= el < pThresholds[i + 1] for el in y] for i in range(len(pThresholds) - 1)])

	for i, mask in enumerate(masks):
		ax.bar(np.array(x)*mask, np.array(y)*mask, color=colors[i], width=1)

	fig.set_size_inches(24, 8)

	#plt.show()
	fig.savefig(outPath)
	plt.close()
	# tell user where the plot was saved
	print >> sys.stderr, outPath

	return 0

### Take a list of indices, a list of scores and a tuple of cutoffs and print a Manhattan plot to outPath
### this plotting function IS threadsafe AFAICT
def manhattanPlot(x, y, outPath, pThresholds=(0.05, 0.01, 0.001), keyID=None):

	xNames = x
	x = range(len(x)) # so you can have a non-numeric x-axis, but they still get plotted in order
	plt.style.use('dark_background')
	fig, ax = plt.subplots()
	y = [-np.log10(p) for p in y] # log-transform p-values
	xlimits = (min(x)-0.5, max(x)+0.5)
	ax.set_xlim(xlimits)
	ax.set_xticks(x, minor=True)
	# hide the default tick labels; this is stupid
	ax.set_xticklabels(x, visible=False)
	# show the correct tick labels
	ax.set_xticklabels(xNames, minor=True, rotation=90, fontsize=8)

	if keyID:
		ax.set_xlabel("position in " + keyID)
	else:
		ax.set_xlabel("predictor column")
	ax.set_ylabel("-log10(P)")

	# add upper and lower thresholds so those bars don't disappear!
	pThresholds = (1,) + pThresholds + (1e-100,)
	pThresholds = [-np.log10(p) for p in pThresholds] # log-transform alpha-values

	# set the threshold colors
	#colors = ["black", "blue", "red", "violet"] # boring colors
	#colors = ["black", "5CBF9F", "1B56CC", "A0242A"] # mpl doesn't take hex!
	#colors = ["black", (0.3608, 0.7490, 0.6235), (0.1059, 0.3373, 0.8000), (0.6275, 0.1412, 0.1647)]# DEEPC colors!
	colors = ["white", (0,1,0.05), (0,1,0.05)]
	# and draw threshold lines
	ax.hlines(pThresholds[:-1], xlimits[0], xlimits[1], colors=colors)

	# for lin
	#masks = [[pThresholds[i] >= el > pThresholds[i + 1] for el in y] for i in range(len(pThresholds) - 1)]
	# for log
	masks = np.array([[pThresholds[i] <= el < pThresholds[i + 1] for el in y] for i in range(len(pThresholds) - 1)])

	for i, mask in enumerate(masks):
		ax.bar(np.array(x)*mask, np.array(y)*mask, color=colors[i], width=1)

	fig.set_size_inches(24, 8)

	#plt.show()
	fig.savefig(outPath)
	plt.close()
	# tell user where the plot was saved
	print >> sys.stderr, outPath

	return 0

### unpack ##all coefficients and their stderrs from a RegressionResults.fit() object
### return a string-indexed pd.DataFrame()
def unpackParams(paramNames, regress, design):
	# get the summary in csv format with lines separated like a file stream
	summary = regress.summary(xname=paramNames).as_csv().split('\n')
	# parse the summary
	reader = list(csv.reader(summary, delimiter=','))
	# strip all the whitespace out of the list of lists
	reader = [[el.strip() for el in line] for line in reader]
	startLine = 10 # line of summary where the data header is
	# cut out the data and headers as a list of lists
	reader = np.array(reader[startLine:startLine+len(paramNames)+1])
	# convert to DataFrame and make the headers indices
	paramDf = pd.DataFrame(reader[1:,1:], index=reader[1:,0], columns=reader[0,1:])
	# add in the nobs
	paramDf["nobs"] = np.sum(design, 0)

	return paramDf


### return t-PDFs and CDFs for non-constant parameters in the passed DataFrame
def tDists(summaryDf):

	# set the distribution domain to a generous 6-sigma using the largest stderr
	xs = np.linspace(np.nanmin(summaryDf["coef"])-6*np.nanmax(summaryDf["std dev"]), np.nanmax(summaryDf["coef"])+6*np.nanmax(summaryDf["std dev"]), num=1000) # make the resolution dynamic as well

	pdfs = dict()  # init dict of PDF lists
	cdfs = dict()  # init dict of PDF lists
	for aa, row in summaryDf.iterrows():
		pdfs[aa] = st.t.pdf(xs, df=row["nobs"] - 1, loc=row["coef"], scale=row["std dev"])
		cdfs[aa] = st.t.pdf(xs, df=row["nobs"] - 1, loc=row["coef"], scale=row["std dev"])

	# return tuple with x-values and dict of distributions keyed on param name
	return (xs, pdfs, cdfs)

### draw superimposed t-distributions for the estimated slopes
### OO API (ax.) is used for thread-safety
def drawDists(xs, dists, depName, colnum, outDir):

	plt.style.use('dark_background')
	fig, ax = plt.subplots()

	for param, dist in dists.items():
		# print >> sys.stderr, dist #TEST
		ax.plot(xs, dist, label=param)

	fig.legend(loc='upper left')
	ax.set_xlabel("{} shift from phylogenetic mean".format(depName))
	ax.set_ylabel("frequency")

	fig.savefig(outDir + "t-PDFs_" + depName + "_col" + str(colnum + 1) + ".pdf")
	plt.close()

	return 0

def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument("-p", "--predictor", type=str, help="predictor values: multiple sequence alignment or trait table filename", required=True)
	parser.add_argument("-f", "--format", default="fasta", help="predictor file format [fasta]; pass 'table' for quantitative mode")
	parser.add_argument('-r', '--response', type=str, help="response values: trait table filename", required=True)
	parser.add_argument("-t", "--tree", type=str, help="tree filename (Newick format)", required=True)
	parser.add_argument("-l", "--lamb-pagel", type=float, default=1, help="Pagel's Lambda [1]") # dammit Pagel, "lambda" is a reserved word!
	parser.add_argument("-s", "--sub_weight", type=float, default=1, help="Substitution rate weight scalar [1]")
	parser.add_argument("-sm", "--sub_matrix", type=str, default="BLOSUM62", help="Pass in a custom substitution rate matrix [BLOSUM62]")
	parser.add_argument("-k", "--key_seq", type=str, help="Name of key sequence on which to index the output columns [None]")
	parser.add_argument("--cpu", type=int, default=cpu_count(), help="Thread count (max # CPUs to use) [{}]".format(cpu_count()))
	parser.add_argument("-b", "--bell_curves", type=float, default=0, help="p-value cutoff below which t-PDF bell curves will be drawn [0]")
	parser.add_argument("-m", "--manhattan", action="store_true", help="Save Manhattan plots with default thresholds [0.05,0.01,0.001]")
	parser.add_argument("-mt", "--manhattan_thresholds", type=tuple, default=(0.10,0.05), help="List of custom thresholds for Manhattan plots")

	global args
	args = parser.parse_args(argv)

	# set mode switch
	quanPredictor = (args.format == "table")

	# Load the continuous response values into a DataFrame in order
	responseTable = makeNiceTraitTable(args.response)

	# tell user what just happened
	print >> sys.stderr, "# MESSAGE: Loaded {} dependent variables for {} taxa".format(responseTable.shape[1], responseTable.shape[0])

	# Load the tree
	tree = Tree(args.tree)
	tree.convert_to_ultrametric(tree_length=1) # IMPORTANT: normalize root-to-leaf distance to 1
	# TODO: look at ways to calculate C from an additive tree. Brownian model may be compromised.

	# Transform branch lengths into a variance-covariance matrix
	# Use the taxa and order provided in responseTable
	# IMPORTANT: orderedLeaves dictates the order of taxa in all matrices and vectors. In this case it should theoretically = responseTable.index
	# IFF they both contain all the same taxa
	phyCovMatrix, orderedLeaves = tree2covMatrix(tree, taxa=responseTable.index ,pagel=args.lamb_pagel) # This stays the same across all columns

	# remake the response table using only taxa that are in the tree
	responseTable = makeNiceTraitTable(args.response, orderedLeaves)
	# and then convert that dataframe to a list of tuples for faster iteration, omitting the index
	responseList = np.array([row for row in responseTable.itertuples(index=False, name=None)])

	# tell user what just happened
	print >> sys.stderr, "# MESSAGE: Generated phylogenetic covariance matrix for {} taxa".format(len(orderedLeaves))

	# make sure there are dirs ready to receive the requested plots
	if args.bell_curves > 0:
		bellDir = "tdists/"
		try:
			os.mkdir(bellDir)
		except:
			pass

	if args.manhattan:
		manhattanDir = "manhattan/"
		try:
			os.mkdir(manhattanDir)
		except:
			pass

	if quanPredictor: # if in quantitative-predictor mode

		# Load the continuous predictor values into a DataFrame and sort it to match the order of the cov matrix
		# this drops any taxa in the table that are missing from the tree
		predictorTable = makeNiceTraitTable(args.predictor, orderedLeaves)
		# and then convert the dataframe to a list of columns for faster iteration, omitting the index
		predictorList = np.array([col for col in predictorTable.T.itertuples(index=False, name=None)])

		# Find predictors that are not in tree or response table, issue warnings
		missingTaxa = set()
		for pTaxon in predictorTable.index:
			if pTaxon not in orderedLeaves:
				#print >> sys.stderr, "# WARNING: {} is missing from the tree and will be dropped from analysis".format(pTaxon) # this is already handled in tree2covMatrix()
				missingTaxa.add(pTaxon)
			if pTaxon not in responseTable.index:
				print >> sys.stderr, "# WARNING: {} is missing from the response table and will be dropped from analysis".format(pTaxon)
				missingTaxa.add(pTaxon)

		def PGLSUnpack(pargs):
			testsResults = list()
			# for each dependent variable
			for depIndex, dependent in enumerate(pargs[0].T):
				# recreate the argument list using only that dependent variable
				dargs = [dependent] + pargs[1:]
				# regress against that variable
				results = pgls(*dargs).fit()
				# in present form I only accept one predictor, so the F-pval and AIC are "all that matter"
				testsResults.append((results.f_pvalue, results.aic))

			pbar.update(1)  # update progress bar
			return testsResults

		# init progress bar
		print >> sys.stderr, "# MESSAGE: Analyzing columns..."
		pbar = tqdm(total=len(predictorList))

		#print responseList
		#print responseList.shape[0]
		#print [predictor for predictor in predictorList]
		#print len(predictorList[0])

		# Processes are split up by column: each process executes 1 GLS regression
		tpool = pool.ThreadPool(args.cpu)
		sumStats = list(tpool.map(PGLSUnpack, [(responseList, predictor, phyCovMatrix) for predictor in predictorList]))
		tpool.close()
		tpool.join()
		pbar.close()

		# prepare the table of p-values for output
		# rows are predictors (genes or traits), columns are response variables
		# put tuples of stats a DataFrame(); index is OG/gene/something else
		outTable1 = pd.DataFrame(sumStats)
		# split the stats into their own columns
		outTable = pd.DataFrame()
		for col in outTable1.columns:
			outTable = pd.concat([outTable,pd.DataFrame(outTable1[col].values.tolist())], axis=1)

		# attach the index IDs and set them to index
		outTable = pd.concat([outTable, pd.DataFrame({"Predictor": predictorTable.T.index})], axis=1).set_index("Predictor")

		# specify stats being reported for each dependent var
		statsNames = ("F-p_value", "F-AIC")
		# iterate over them to get the column names in order (to match the DataFrame() generated above)
		outTable.columns = [depName + '.' + statName for depName in responseTable.columns for statName in statsNames]

	else: # if in CPGLS mode

		# Load in AA alignment
		alignment = AlignIO.read(args.predictor, format=args.format)

		# Load the substitution matrix and normalize if necessary
		if args.sub_matrix.find("BLOSUM") > -1: # if it's a BLOSUM matrix
			try: # look for it in the blosum module
				subMatrix = blosum.submatrix(int(''.join([i for i in args.sub_matrix if i.isdigit()])))
			except:
				print >> sys.stderr, "# ERROR: {} is not part of the blosum submatrix module!".format(args.sub_matrix)
				sys.exit(1)
		else: # load it from file
			subMatrix = np.array(pd.read_table(args.sub_matrix, sep='\t'))

		if not np.allclose(subMatrix, subMatrix.T): # NOTE that the 'is' operator won't work with np.array()!
			print >> sys.stderr, "# ERROR: Substitution rate matrix is asymmetric, exiting"
			sys.exit(1)
		subMatrixNorm = normalizeRowWise(subMatrix, 1)
		if not np.allclose(subMatrix, subMatrixNorm):
			print >> sys.stderr, "# WARNING: Substitution rate matrix was not properly normalized, but now it is"
			subMatrix = subMatrixNorm

		# Find seqs that are not in tree or response table, issue warnings
		missingTaxa = set()
		for record in alignment:
			if record.id not in orderedLeaves:
				print >> sys.stderr, "# WARNING: {} is missing from the tree and will be dropped from analysis".format(record.id)
				missingTaxa.add(record.id)
			if record.id not in responseTable.index:
				print >> sys.stderr, "# WARNING: {} is missing from the response table and will be dropped from analysis".format(record.id)
				missingTaxa.add(record.id)

		# Drop those missing taxa
		keepRecords = dict()
		for record in alignment:
			if record.id not in missingTaxa:
				keepRecords[record.id] = record
		keepRecordsList = [keepRecords[taxon] for taxon in orderedLeaves] # get the seqs in order!

		# Get the root phenotypes, aka phylomeans, for fixing the intercept
		# unpack the response vars into taxon-keyed dicts
		responseDicts = [dependent[1].to_dict() for dependent in responseTable.T.iterrows()]
		# get the phylomeans and store them to a list of intercepts
		# the below pulls out the _last_ value for a tree traversed in post-order
		intercepts = [ancR(tree, response)[-1] for response in responseDicts]

		# store the culled alignment as an array, sequences row-wise
		alignArray = np.array([list(record) for record in keepRecordsList], np.character)
		# then split in into column vectors
		colVectors = np.hsplit(alignArray, alignArray.shape[1])

		if args.key_seq:  # map the PP scores to key sequence positions
			# get key columns of the alignment
			keyCols = key_alignment_columns(alignment, args.key_seq)
			colVectors = [colVectors[i-1] for i in keyCols]

		### Multiple-arg unpacker for the ThreadPool. Dependent variable names and bell curve directory are passed by default
		def CPGLSUnpack(pargs, depNames=responseTable.columns):
			# init master list for all significance test results across dependents
			testsResults = list()
			# for each dependent variable
			for depIndex, dependent in enumerate(pargs[0].T):
				# get the name of the dependent variable
				depName = depNames[depIndex]
				intercept = intercepts[depIndex]
				# recreate the argument list using only that dependent variable
				dargs = [dependent, intercept] + pargs[2:]
				#(responseList, intercepts, phyCovMatrix, subMatrix, column, args.sub_weight, colnum)
				#(responseList, intercept, phyCovMatrix, subMatrix, colVector, subWeight=1, colnum=0)
				# regress against that variable
				regDesign = cpgls(*dargs)

				# can get the whitened data for further manipulation
				#print "whitened endog"
				#print pd.DataFrame(regDesign.wendog).to_csv(sep='\t') #TEST, print whitened data
				#print "whitened exog"
				#print pd.DataFrame(regDesign[0].wexog).to_csv(sep='\t')  # TEST, print whitened data

				results = regDesign.fit()  # fit the regression, SILENTLY
				# if we want to see the summary of every regression
				#print results.summary()

				#print "FP " + str(results.f_pvalue) #TEST

				if ~np.isnan(results.f_pvalue): # if not a conserved column

					# do pairwise 2-sample t-tests and get the p-values
					multicorr = 'hs' # the multiple-test correction to be used
					pairwiseT = results.t_test_pairwise(term_name='aa', method=multicorr).result_frame
					pairTDict = pairwiseT["pvalue-"+multicorr].T.to_dict()
					min_pairT = np.nanmin(pairTDict.values())
					#pairTDict = str(pairTDict) # safer for .to_csv()?

					# Do a within/among comparison of variances
					anovaF_pvalue = sm.stats.anova_lm(results, typ=2)["PR(>F)"]["aa"]

					# Do normality tests for diagnostic purposes
					# omni:
					omniNorm_pvalue = sms.omni_normtest(results.resid)[1]
					# Jarque-Bera:
					#jbNorm_pvalue = sms.jarque_bera(results.resid)[1]

					# collate all the p-values, including the F-pVal for the linear model
					pDictList = [results.f_pvalue, anovaF_pvalue, min_pairT, pairTDict, omniNorm_pvalue]

				else: pDictList = [None] * 5 # return a null list; must be same length as pDictList above!

				# append the set of p-values to the list of lists for different response variables
				testsResults = testsResults + pDictList

				# if the column makes the specified p-cutoff to draw bell curves
				'''if pDictList[2] and pDictList[2] <= args.bell_curves:
					# generate discrete series for the t-distributions
					xs, tPDFs, tCDFs = tDists(summary)
					# draw t-distributions identified by dependent and column
					drawDists(xs, tPDFs, depName, pargs[5], bellDir)'''

			pbar.update(1) # update progress bar
			# return minimum p-val and p-val dict for pool map, all should be Bonferroni-corrected beforehand
			return testsResults

		# init progress bar#
		print >> sys.stderr, "# MESSAGE: Analyzing columns..."
		pbar = tqdm(total=len(colVectors))

		#for i in range(len(colVectors)): #TEST for error handling in the ThreadPool
		#	print >> sys.stderr, CPGLSUnpack([responseList, intercepts, phyCovMatrix, subMatrix, colVectors[i], args.sub_weight, i])

		with warnings.catch_warnings():
			warnings.simplefilter("ignore")
			# Processes are split up by column: each process executes 1 GLS regression
			tpool = pool.ThreadPool(args.cpu)
			sumStats = list(tpool.map(CPGLSUnpack, [[responseList, intercepts, phyCovMatrix, subMatrix, column, args.sub_weight, colnum] for colnum, column in enumerate(colVectors)])) # split alignment array into columns and feed them to pgls()
			tpool.close()
			tpool.join()
			pbar.close()

		# prepare the table of p-values for output
		outTable = pd.DataFrame({"Sites": [i+1 for i in range(len(colVectors))]}).set_index("Sites") # init table with a site column
		# specify stats being reported for each dependent var
		statsNames = ("model_F-p", "anova_F-p", "min_t-p_HolmSidak", "pairwise_t-p_HolmSidak", "omni_X2-p") # to return the min p-val and all the p-vals
		# create data columns in the same order as the lists returned by the ThreadPool
		statsColumns = [depName + '.' + statName for depName in responseTable.columns for statName in statsNames]
		for i, statCol in enumerate(statsColumns):
			outTable[statCol] = [el[i] for el in sumStats]

	if args.manhattan:
		# NOTE that this will generally crash any parent script calling this main() function,
		# because matplotlib fails to hand plotting resources back to the parent process.

		# get all numeric columns from outTable
		numerics = ['int16', 'int32', 'int64', 'float16', 'float32', 'float64']
		outNumerics = outTable.select_dtypes(include=numerics)

		# draw plots in parallel
		print >> sys.stderr, "# MESSAGE: Drawing Manhattan plots..."
		def manhattan_unpack(pargs):
			manhattanPlot(*pargs)
			return 0

		#print >> sys.stderr, manhattan_unpack(*[(outTable.index, list(outNumerics.T.ix[0]), manhattanDir+str(outNumerics.T.index[0])+".pdf", args.manhattan_thresholds, args.key_seq)]) #TEST

		tpool = pool.ThreadPool(args.cpu)
		#tpool = pool.ThreadPool(1)
		# list comp plots each numeric column against the index and names the files for those numeric columns
		# with itertuples(), element 0 of the row is the index name and [1:] are the actual values
		tpool.map(manhattan_unpack, [(outTable.index, list(y[1:]), manhattanDir+str(y[0])+".pdf", args.manhattan_thresholds, args.key_seq) for y in outNumerics.T.itertuples(name=None)])
		tpool.close()
		tpool.join()

	# return the results DataFrame to any calling process
	return outTable


if __name__ == "__main__":
	# when called from command line, print results to stdout
	print main(sys.argv[1:], sys.stdout).to_csv(sep = '\t')