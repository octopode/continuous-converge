#!/usr/bin/python

import pcoc_cont_scenarios as target

import os


'''
def test_init_tree():

    # input Newick data
    newick = "((((A:2.0,B:2.0):1.0,C:3.0):1.0,((D:2.0,E:2.0):1.0,F:3.0):1.0):4.0,(G:4.0,H:4.0):4.0);"

    t = target.init_tree(newick)
    print t
    print len([n for n in t.traverse("postorder")])

    assert False
'''

# test ancestral reconstruction function
def test_ancR():

    # input ete3 tree and trait dict
    newick = "((((A:2.0,B:2.0):1.0,C:3.0):1.0,((D:2.0,E:2.0):1.0,F:3.0):1.0):4.0,(G:4.0,H:4.0):4.0);"
    tree = target.init_tree(newick)
    tipTraits = {"A":    1,
                 "B":    2,
                 "C":    3,
                 "D":    4,
                 "E":    5,
                 "F":    6,
                 "G":    7,
                 "H":    8}

    # output node trait list
    output = [1, 2, 1.5, 3, 1.875, 4, 5, 4.5, 6, 4.875, 3.375, 7, 8, 7.5, 5.4375]

    # diagnostic output
    # print [node.name for node in tree.traverse("postorder")]
    # target.traitTree(tree, output, target.wavelength_to_rgb, "/data", float=3)

    # uncomment to print result
    #print target.ancR(tree, tipTraits); nodeTraits = False

    logger = target.logging.getLogger("test")
    assert target.ancR(tree, tipTraits, logger) == output

def test_binBy_nBins():

    input = [1,2,3,4,5,6,7,8,9,10,11,12]
    numBins = 3

    output = ([4.5, 8.5], [2.5, 6.5, 10.5])

    logger = target.logging.getLogger("test")
    assert target.binBy(input, logger, fixNumBins = numBins) == output

def test_binBy_binWidth():
    input = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    binWidth = 4

    output = ([4.5, 8.5], [2.5, 6.5, 10.5])

    logger = target.logging.getLogger("test")
    assert target.binBy(input, logger, fixBinWidth = binWidth) == output

def test_discretize():

    # input trait values
    traits = [1, 2, 1.5, 3, 1.875, 4, 5, 4.5, 6, 4.875, 3.375, 7, 8, 7.5, 5.4375]
    cutoffs = (1.25, 1.6875, 1.9375, 2.5, 3.1875, 3.6875, 4.25, 4.6875, 4.9375, 5.21875, 5.71875, 6.5, 7.25, 7.75)

    # output discretized trait table
    discrete = [[False,False,False,False,False,False,False,False,False,False,False,False,False,False],
                [True,True,True,False,False,False,False,False,False,False,False,False,False,False],
                [True,False,False,False,False,False,False,False,False,False,False,False,False,False],
                [True,True,True,True,False,False,False,False,False,False,False,False,False,False],
                [True,True,False,False,False,False,False,False,False,False,False,False,False,False],
                [True,True,True,True,True,True,False,False,False,False,False,False,False,False],
                [True,True,True,True,True,True,True,True,True,False,False,False,False,False],
                [True,True,True,True,True,True,True,False,False,False,False,False,False,False],
                [True,True,True,True,True,True,True,True,True,True,True,False,False,False],
                [True,True,True,True,True,True,True,True,False,False,False,False,False,False],
                [True,True,True,True,True,False,False,False,False,False,False,False,False,False],
                [True,True,True,True,True,True,True,True,True,True,True,True,False,False],
                [True,True,True,True,True,True,True,True,True,True,True,True,True,True],
                [True,True,True,True,True,True,True,True,True,True,True,True,True,False],
                [True,True,True,True,True,True,True,True,True,True,False,False,False,False]]
    output = target.pd.DataFrame(columns = cutoffs, data = discrete)

    assert target.discretize(traits, cutoffs).equals(output)

def test_uniqueScenarios():

    # input discrete trait table. Columns 3, 4, 5 are duplicate of 0, 1, 2.
    cutoffs = (1,   2,  3,  4,  5,  6)
    discrete = [[False, False,  False,  False,  False,  False],
                [False, False,  True,   False,  False,  True],
                [False, True,   True,   False,  True,   True],
                [False, False,  False,  False,  False,  False],
                [False, False,  True,   False,  False,  True],
                [False, True,   True,   False,  True,   True]]
    input = target.pd.DataFrame(columns=cutoffs, data=discrete)

    # output consolidated dataframe
    avgCutoffs = (2.5,  3.5,    4.5)
    conDiscrete = [[False,False,False],
                    [False,False,True],
                    [False,True,True],
                    [False,False,False],
                    [False,False,True],
                    [False,True,True]]
    output = target.pd.DataFrame(columns=avgCutoffs, data=conDiscrete)

    logger = target.logging.getLogger("test")
    assert target.uniqueScenarios(input, logger).equals(output)

def test_convergentRootFilter():

    # input Boolean dataframe
    # 2/3 scenarios here have a convergent root
    cutoffs = (2.5,  3.5,    4.5)
    discrete = [[False, False, False],
                   [False, False, True],
                   [False, True, True],
                   [False, False, False],
                   [False, False, True],
                   [False, True, True]]
    input = target.pd.DataFrame(columns=cutoffs, data=discrete)

    # output consolidated dataframe
    filteredCutoffs = (2.5,)
    filteredDiscrete = [[False],
                [False],
                [False],
                [False],
                [False],
                [False]]
    output = target.pd.DataFrame(columns=filteredCutoffs, data=filteredDiscrete)

    logger = target.logging.getLogger("test")
    assert target.convergentRootFilter(input, logger).equals(output)

def test_scenarioStrings():

    # input ete3 tree and trait dict
    newick = "((((A:2.0,B:2.0):1.0,C:3.0):1.0,((D:2.0,E:2.0):1.0,F:3.0):1.0):4.0,(G:4.0,H:4.0):4.0);"
    tree = target.init_tree(newick)

    # input discretized trait table
    cutoffs = (1.25, 1.6875, 1.9375, 2.5, 3.1875, 3.6875, 4.25, 4.6875, 4.9375, 5.21875, 5.71875, 6.5, 7.25, 7.75)
    discrete = [[False, False, False, False, False, False, False, False, False, False, False, False, False, False],
                [True, True, True, False, False, False, False, False, False, False, False, False, False, False],
                [True, False, False, False, False, False, False, False, False, False, False, False, False, False],
                [True, True, True, True, False, False, False, False, False, False, False, False, False, False],
                [True, True, False, False, False, False, False, False, False, False, False, False, False, False],
                [True, True, True, True, True, True, False, False, False, False, False, False, False, False],
                [True, True, True, True, True, True, True, True, True, False, False, False, False, False],
                [True, True, True, True, True, True, True, False, False, False, False, False, False, False],
                [True, True, True, True, True, True, True, True, True, True, True, False, False, False],
                [True, True, True, True, True, True, True, True, False, False, False, False, False, False],
                [True, True, True, True, True, False, False, False, False, False, False, False, False, False],
                [True, True, True, True, True, True, True, True, True, True, True, True, False, False],
                [True, True, True, True, True, True, True, True, True, True, True, True, True, True],
                [True, True, True, True, True, True, True, True, True, True, True, True, True, False],
                [True, True, True, True, True, True, True, True, True, True, False, False, False, False]]
    input = target.pd.DataFrame(columns=cutoffs, data=discrete)

    output = {1.25: '14,13,12,11,10,9,8,7,6,5,4,3,2,1',
              1.6875: '14,13,12,11,10,9,8,7,6,5,4,3/1',
              4.9375: '14,13,12,11/8/6',
              7.25: '13,12',
              3.6875: '14,13,12,11/9,8,7,6,5',
              6.5: '13,12,11',
              7.75: '12',
              4.6875: '14,13,12,11/9,8/6',
              5.71875: '13,12,11/8',
              5.21875: '14,13,12,11/8',
              4.25: '14,13,12,11/9,8,7,6',
              1.9375: '14,13,12,11,10,9,8,7,6,5/3/1',
              2.5: '14,13,12,11,10,9,8,7,6,5/3',
              3.1875: '14,13,12,11,10,9,8,7,6,5'}

    assert target.scenarioStrings(input, tree) == output

# at this point, reconnect the filters and run the program
# then, refactor the PCOC scripts.

# test the site-determine module
def disabled_test_det_main():

    # directory this test_...py script is in
    dir_path = os.path.dirname(os.path.realpath(__file__))

    treeFile = dir_path + "/test_files/ABCDEFGH.tree"
    aliFile = dir_path + "/test_files/Converge7state_realAcids_all.fasta"
    traitFile = dir_path + "/test_files/IntegerDepths_8taxa.tsv"
    testScenario = "14,13,12,11/9,8/6"

    outDir = dir_path + "/test_files/test_output"

    detArgs = target.argparse.Namespace(
                    # dynamic values
                    manual_mode = testScenario,
                    # static values
                    tree = treeFile,
                    ali = aliFile,
                    output_dir = outDir,
                    # default values
                    ph = None,
                    debug = False,
                    LD_LIB = "",
                    CATX_est = 10,
                    gamma = False,
                    inv_gamma = False,
                    max_gap_allowed = 5,
                    max_gap_allowed_in_conv_leaves = 5,
                    no_cleanup = False,
                    filter_t = .99,
                    filter_t_pcoc = -1,
                    filter_t_pc = -1,
                    filter_t_oc = -1,
                    plot = False,
                    reorder = False,
                    svg = False,
                    no_cleanup_fasta = False
                )

    target.det.main(args = detArgs)

    assert False

def disabled_test_det_mk_detect():
    # CAUTION: this one takes awhile

    # input
    manual_mode_nodes = {"T":[14, 9, 6],
                         "C":[13, 12, 11, 8]}

    # directory this test_...py script is in
    dir_path = os.path.dirname(os.path.realpath(__file__))

    tree_filename = dir_path + "/test_files/ABCDEFGH.tree"
    # try loading the Newick into RAM; note that it is no longer actually a filename
    tree_filename = target.Tree(tree_filename).write()
    ali_filename = dir_path + "/test_files/Converge7state_realAcids_all.fasta"
    OutDirName = dir_path + "/test_files/test_output"

    tempDirs = target.det.get_temp_dirs(ali_filename, OutDirName)

    ali = target.AlignIO.read(ali_filename, "fasta")

    # output
    outDict = {'PCOC': {0: 0.98776456995135431, 1: 0.770431306451275, 2: 0.89761946663107584, 3: 0.0991179423639134,
                        4: 0.11497691785505021, 5: 0.65316350744083285, 6: 0.54772777783491067, 7: 0.56506668895596113,
                        8: 0.082620720326110342, 9: 0.00058786069163569701},
               'OC': {0: 0.48526744245668019, 1: 0.44903619735509348, 2: 0.47398015275966554, 3: 0.47704571896670039,
                      4: 0.47590655215689426, 5: 0.50825351823505061, 6: 0.50484097383275128, 7: 0.50887462012502926,
                      8: 0.50159629752294632, 9: 0.48475219546269044},
               'Sites': {0: 1, 1: 2, 2: 3, 3: 4, 4: 5, 5: 6, 6: 7, 7: 8, 8: 9, 9: 10},
               'PC': {0: 0.99766551234963974, 1: 0.90931033001525507, 2: 0.96521065373097248, 3: 0.16643203866962952,
                      4: 0.18284000859464614, 5: 0.64098650008475289, 6: 0.54361244613942983, 7: 0.55684601648332299,
                      8: 0.0890370485404164, 9: 0.000727989717059959},
               'Indel_prop': {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0, 6: 0.0, 7: 0.0, 8: 0.0, 9: 0.0},
               'Indel_prop(ConvLeaves)': {0: 0.0, 1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0, 5: 0.0, 6: 0.0, 7: 0.0, 8: 0.0,
                                          9: 0.0}}
    # it's important to sort the columns or DataFrame.equals() won't work
    output = target.pd.DataFrame.from_dict(outDict).sort_index(axis=1)

    logger = target.logging.getLogger("test")
    assert target.det.mk_detect(manual_mode_nodes, tree_filename, ali, logger, *tempDirs).sort_index(axis=1).equals(output)

def test_trawlErrorCurves():

