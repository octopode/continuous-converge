import sys
from ete3 import Tree, NodeStyle, TreeStyle, TextFace, CircleFace

### Number tree nodes for consistent reference
def init_tree(nf):
    t = Tree(nf)

    for i, n in enumerate(t.traverse("postorder")):
        n.add_features(ND = i)

    return t

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
#def testTrees(tree, scenarios, outDir, treePath, floatSwitch=0):
def testTrees(scenarios, outDir, treePath, floatSwitch=0):

    ### Draw test trees. This is a modified version of the test routine in pcoc_num_tree.py, stuffed in a for loop
    for cutoff in sorted(scenarios.keys()):
        tree = init_tree(treePath)
        # not mucking with additive trees yet; ultrametrize the tree and normalize to length 1
        #tree.convert_to_ultrametric(tree_length=1)
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