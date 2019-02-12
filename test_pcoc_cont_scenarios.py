#!/usr/bin/python

import pcoc_cont_scenarios as target


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

    assert target.ancR(tree, tipTraits) == output

def test_binBy_nBins():

    input = [1,2,3,4,5,6,7,8,9,10,11,12]
    numBins = 3

    output = ([4.5, 8.5], [2.5, 6.5, 10.5])

    assert target.binBy(input, fixNumBins = numBins) == output

def test_binBy_binWidth():
    input = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    binWidth = 4

    output = ([4.5, 8.5], [2.5, 6.5, 10.5])

    assert target.binBy(input, fixBinWidth = binWidth) == output

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

    assert target.uniqueScenarios(input).equals(output)

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

    assert target.convergentRootFilter(input).equals(output)

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