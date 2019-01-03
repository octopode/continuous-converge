#!/usr/bin/env python
#
# pdb_encode_sitewise.py v1 2018-12-06

'''pdb_site_identity.py  last modified 2018-09-12

pdb_site_identity.py -a mox_all.aln -s DOPO_HUMAN -p 4zel.pdb > 4zel_w_scores.pdb

within a PDB file, fields for atoms are:
Record name	Residue	Position as X Y Z	Atom serial number	Occupancy	Atom	Chain	Residue number
such as:
ATOM      1  N   PRO A  46       0.739  40.031  44.896  1.00 0.842           N

here, the temperatureFactor will be replaced with the identity %,
where bins are:
0, 50, 60, 70, 80, 90, 95, 98, 100%

proteins or ligands not in the alignment are coded as -1, and colored "null"

residue number must match the alignment, not the position in the model
meaning even if the PDB file starts with residue 20, if the first 19 were
disordered or cleaved, the sequence still must start with residue 1
'''

import sys
import time
import argparse
from collections import Counter,defaultdict
from numpy import log, mean
import pandas as pd
from Bio import AlignIO


'''def tables2dict(csvpaths):
	# Take a list of paths to CSV files and return a dict of DataFrames keyed on the extensionless basename of the file
	data = dict()
	for fname in csvpaths:
		keyname = os.path.splitext(os.path.basename(fname))[0]
		try:
			data[keyname] = pd.read_csv(fname)
		except:
			print >> sys.stderr, "WARNING: DATA FILE {} NOT FOUND".format(fname)

	return data'''


def key_alignment_columns(alignment, key_seqid):
	'''return list of columns of an alignment that are not gaps in the key sequence'''

	al_length = alignment.get_alignment_length()

	# Get the key sequence
	for seqrec in alignment:
		if seqrec.id == key_seqid:
			keyseq = seqrec.seq

	mappedCols = [None] * len(str(keyseq).replace('-',''))
	keyindex = 0
	for i in range(al_length):  # only iterate over columns that are not gaps in target seq
		if keyseq[i] != "-":  # meaning anything except gaps
			mappedCols[keyindex] = i+1 # map the alignment column to the key column. Gotta shift the index!
			keyindex += 1

	print >> sys.stderr, mappedCols
	return mappedCols


def check_occupancy(alignment, key_seqid):
	'''Utility function to assess whether the passed reference seqid is a good choice.
	In other words, do gaps in the reference have high occupancy in the alignment, or are they rare insertions?'''

	al_length = alignment.get_alignment_length()
	num_taxa = len(alignment)

	# Get the key sequence
	for seqrec in alignment:
		if seqrec.id == key_seqid:
			keyseq = seqrec.seq

	gapOccupancy = [None] * keyseq.count('-')
	gapindex = 0
	for i in range(al_length):  # here, iterate over gaps in key only!
		if keyseq[i] == "-":  # meaning only gaps!
			occupancy = 0 # init
			for seqrec in alignment: # for each taxon
				if seqrec.seq[i] != '-': # if that position is occupied
					occupancy += 1 # count it!
			occupancy /= num_taxa # and then normalize occupancy by the number of taxa
			gapOccupancy[gapindex] = occupancy
			gapindex += 1

	print >> sys.stderr, "# Key seq {} has {} gaps with a mean occupancy of {}%".format(key_seqid, len(gapOccupancy), mean(gapOccupancy)*100)

	return gapOccupancy



def get_chains_only(defaultchain, seqidlist, pdbfile):
	'''read PDB file and return a dict where key is chain and value is sequence ID'''
	keepchains = {} # dict where key is chain and value is seqid
	print >> sys.stderr, "# Reading chain from PDB {}".format(pdbfile)
	for line in open(pdbfile,'r'):
		record = line[0:6].strip()
		# get relevant chains that match the sequence, in case of hetero multimers
		if record=="DBREF":
			defaultchain = False
			proteinid = line[42:56].strip()
			for seqid in seqidlist:
				if seqid.find(proteinid)>-1:
					chaintarget = line[12]
					print >> sys.stderr, "### keeping chain {} for sequence {}".format( chaintarget, proteinid )
					keepchains[chaintarget] = proteinid
	if defaultchain: # meaning nothing was found, use default and single sequence
		keepchains[defaultchain] = seqidlist[0]
	return keepchains


def make_output_script(scoredict, keepchains, colorscript, basecolor="red"):
	'''from the identity calculations, print a script for PyMOL'''
	colors = {} # key is colorscheme name, value is list of colors
	colors["red"] = [ [0.63,0.63,0.63] , [0.73,0.55,0.55] , [0.75,0.47,0.47],
				   [0.77,0.38,0.38] , [0.79,0.29,0.29] , [0.82,0.21,0.21],
				   [0.84,0.13,0.13]    , [0.88,0,0]     , [1,0,0.55] ]
	colors["green"] = [ [0.63,0.63,0.63] , [0.50,0.68,0.56] , [0.42,0.71,0.53],
				   [0.35,0.74,0.49] , [0.26,0.77,0.44] , [0.19,0.80,0.41],
				   [0.12,0.83,0.37] , [0.01,0.87,0.31] , [1,0,0.55] ]
	colors["blue"] = [ [0.63,0.63,0.63] , [0.50,0.58,0.68] , [0.42,0.55,0.71],
				   [0.35,0.52,0.73] , [0.28,0.49,0.76] , [0.20,0.46,0.80],
				   [0.12,0.43,0.83] , [0.00,0.38,0.87] , [1,0,0.55] ]
	colors["yellow"] = [ [0.63,0.63,0.63] , [0.66,0.67,0.51] , [0.68,0.71,0.43],
				   [0.70,0.73,0.36] , [0.72,0.76,0.28] , [0.74,0.80,0.20],
				   [0.76,0.83,0.12] , [0.79,0.87,0.00] , [1,0,0.55] ]
	insuf_color = [0.75, 0.75, 0.58]

	#binvalues = [0.0, 50.0 ,60.0 ,70.0 ,80.0 ,90.0 ,95.0 ,98.0 ,100, 101]
	binvalues = [0.0, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1.0]

	print >> sys.stderr, "# Generating PyMOL script {}".format(colorscript)
	# begin printing commands for PyMOL script
	with colorscript as cs:
		print >> cs, "hide everything"
		#print >> wayout, "bg white"
		print >> cs, "show cartoon"
		print >> cs, "set_color colordefault, [{}]".format( ",".join(map(str,insuf_color)) )
		print >> cs, "color colordefault, all"
		# make commands for each color
		for color, colorlist in colors.iteritems():
			for i,rgb in enumerate(colorlist):
				colorname = "{}{}".format( color, binvalues[i] )
				print >> cs, "set_color {}, [{}]".format( colorname, ",".join(map(str,rgb)) )

		# make commands for each chain
		for chain in keepchains.iterkeys(): # keys are chain letters, values are seq IDs
			pctgroups = defaultdict(list) # key is percent group, value is list of residues
			# for each residue, assign to a bin
			for residue in scoredict[keepchains[chain]].iterkeys():
				residuescore = scoredict[keepchains[chain]].get(residue,0.00)
				for i,value in enumerate(binvalues[:-1]):
					upper = binvalues[i+1]
					if residuescore < upper:
						pctgroups[value].append(residue)
						break
				# should not need an else
			# assign whole chain to lowest color, then build up
			print >> cs, "color {}{}, chain {}".format( basecolor, binvalues[0], chain )
			# for each bin, make a command to color all residues of that bin
			for i,value in enumerate(binvalues[:-1]):
				if i==0: # long lists apparently crash the program, so skip
					continue
				binname = "{}score_grp_{}_{}".format( value, i+1, chain )
				resilist = map(str,pctgroups[value])
				binresidues = ",".join(resilist)
				print >> cs, "select {}, (chain {} & resi {})".format( binname, chain, binresidues )
				print >> cs, "color {}{}, {}".format( basecolor, value, binname )
	# no return


def rewrite_pdb(pdbfile, seqidlist, scoredict, wayout, forcerecode=False):
	print >> sys.stderr, "# Reading PDB from {}".format(pdbfile)
	atomcounter = 0
	residuecounter = {}
	keepchains = {} # dict where key is chain and value is seqid
	defaultchain = True # flag for whether DBREF occurs at all
	for line in open(pdbfile,'r'):
		# records include:
		# HEADER TITLE COMPND SOURCE AUTHOR REVDAT JRNL REMARK
		# DBREF SEQRES HET HETNAM FORMUL HELIX SHEET SSBOND LINK CISPEP SITE ATOM CONECT

		#COLUMNS        DATA  TYPE    FIELD        DEFINITION
		#-------------------------------------------------------------------------------------
		# 1 -  6        Record name   "ATOM  "
		# 7 - 11        Integer       serial       Atom  serial number.
		#13 - 16        Atom          name         Atom name.
		#17             Character     altLoc       Alternate location indicator.
		#18 - 20        Residue name  resName      Residue name.
		#22             Character     chainID      Chain identifier.
		#23 - 26        Integer       resSeq       Residue sequence number.
		#27             AChar         iCode        Code for insertion of residues.
		#31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
		#39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
		#47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
		#55 - 60        Real(6.2)     occupancy    Occupancy.
		#61 - 66        Real(6.2)     tempFactor   Temperature  factor.
		#77 - 78        LString(2)    element      Element symbol, right-justified.
		#79 - 80        LString(2)    charge       Charge  on the atom.

		record = line[0:6].strip()
		# get relevant chains that match the sequence, in case of hetero multimers
		if record=="DBREF":
			defaultchain = False
			proteinid = line[42:56].strip()
			for seqid in seqidlist:
				if seqid.find(proteinid)>-1:
					chaintarget = line[12]
					print >> sys.stderr, "### keeping chain {} for sequence {}".format( chaintarget, proteinid )
					keepchains[chaintarget] = proteinid
		# DBREF lines should come before ATOM lines, so for all other lines, check for ATOM or not
		if record=="ATOM": # skip all other records
			chain = line[21]
			residue = int( line[22:26] )
			if defaultchain or forcerecode or chain in keepchains: # default chain means take all, or use chain A
				atomcounter += 1
				if defaultchain or forcerecode: # assume only one seqid
					residuescore = scoredict[seqidlist[0]].get(residue,0.00)
				else:
					residuescore = scoredict[keepchains[chain]].get(residue,0.00)
				if residuescore:
					residuecounter[residue] = True
			else: # meaning in another chain, so color as insufficient
				residuescore = -1
			newline = "{}{:6.2f}{}".format( line[:60], residuescore, line[66:].rstrip() )
			print >> wayout, newline
		else: # this will also print DBREF lines
			print >> wayout, line.strip()
	if atomcounter:
		print >> sys.stderr, "# Recoded values for {} atoms in {} residues".format(atomcounter, len(residuecounter) )
	else:
		print >> sys.stderr, "# NO CHAINS FOUND MATCHING SEQ ID {}, CHECK NAME {}".format( seqid, proteinid )



def main(argv, wayout):
	if not len(argv):
		argv.append('-h')
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument("-d", "--data-table", nargs="*", help="Tab-delimited data to map to structure", required=True)
	parser.add_argument("-a","--alignment", nargs="*", help="multiple sequence alignment", required=True)
	#parser.add_argument("-c","--conservation", action="store_true", help="calculate sitewise positional conservation (CT model)")
	parser.add_argument("-f","--format", default="fasta", help="alignment format [fasta]")
	#parser.add_argument("-g","--gap-cutoff", default=0.5, type=float, help="minimum fraction of non-gap characters per site, else is called unconserved [0.5]")
	parser.add_argument("-p","--pdb", help="PDB format file", required=True)
	parser.add_argument("-s","--sequence", nargs="*", help="key seqid, must match DBREF for PDB", required=True)
	parser.add_argument("-w","--write-script", action="store_true", help="write to script instead of recoding PDB file")
	parser.add_argument("--base-color", default="red", help="default color gradient [red,yellow,green,blue]")
	parser.add_argument("--default-chain", default="A", help="default letter of chain, if DBREF for the sequence cannot be found in PDB [A]")
	parser.add_argument("--force-recode", action="store_true", help="force recoding regardless of chain")
	args = parser.parse_args(argv)

	# Runtime checks:
	if len(args.data_table) != len(args.alignment) != len(args.sequence):
		print >> sys.stderr, "ERROR: {} DATASETS FOR {} ALIGNMENTS FOR {} CHAINS, MUST BE EQUAL, CHECK -a AND -s".format(len(args.data_table), len(args.alignment), len(args.sequence)), time.asctime()

	if len(args.sequence) > len(set(args.sequence)):
		print >> sys.stderr, "ERROR: NON UNIQUE NAMES FOR SEQUENCES, CHECK -s"

	# Load in alignments
	alignments = [AlignIO.read(fname, format=args.format) for fname in args.alignment]

	# Load in a list of dataframes that are scores indexed across alignment sites. 1 dataframe per alignment.
	# Note that dataset-alignment-chain matching is positional. They must be passed to the script in corresponding order like:
	# -d A.tab B.tab C.tab -a A.fa B.fa C.fa -s A B C
	dataframes = [pd.read_csv(fname, sep='\t') for fname in args.data_table]

	if not dataframes: # if the program fails to load *any* sitewise data
		sys.exit("ERROR: NO SITEWISE DATA LOADED, CHECK -d AND YOUR INPUT FILES")
	else: # proceed
		scoredict = defaultdict() # init scoredict
		if args.write_script: # write output PyMOL script with color commands
			for i, algtscores in enumerate(dataframes): # for each chain
				#print >> sys.stderr, algtscores #TEST
				algtscores = algtscores["PCOC"] #TEST - pulls one column out of the score table as a list
				# assess choice of reference
				check_occupancy(alignments[i], args.sequence[i])
				# map scores to the reference seq
				keycolumns = key_alignment_columns(alignments[i], args.sequence[i])
				keyscores = {keycolumns.index(j)+1:algtscores[j-1] for j in keycolumns} # gotta make this a dict to work with WRFs original script function. Gotta shift index back!
				scoredict[args.sequence[i]] = keyscores

			refchains = get_chains_only(args.default_chain, args.sequence, args.pdb)
			make_output_script(scoredict, refchains, sys.stdout, args.base_color) # put sys.stdout where the script fname previously was
			# scoredict key is target seqid, value is dict of scores
		else: # recode PDB beta-factors
			sys.exit("Sorry, PDB recoding is not supported yet")
			#rewrite_pdb(args.pdb, args.sequence, conservedict, wayout, args.force_recode) # pass score dict into file

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
