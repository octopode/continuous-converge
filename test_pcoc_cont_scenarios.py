#!/usr/bin/python

import os
import sys
import pcoc_cont_scenarios as target

# test ancestral reconstruction function
def test_ancR():

    # input ete3 tree and trait dict
    newick = "((C:2,((E:1.5,D:1.5):1.5):1):1,B:1,A:1);"
    tree = target.init_tree(newick)
    tipTraits = {"A":    1,
                 "B":    2,
                 "C":    3,
                 "D":    4,
                 "E":    5}

    print tree
    
    ndDict = {node.ND: node.name for node in tree.traverse("postorder")}

    print ndDict

    # output node trait list
    nodeTraits = []

    assert target.ancR(tree, tipTraits) == nodeTraits
