import sys
from ete3 import Tree, NodeStyle, TreeStyle, TextFace, CircleFace, RectFace

### Number tree nodes for consistent reference
def init_tree(nf):
    t = Tree(nf)

    for i, n in enumerate(t.traverse("postorder")):
        n.add_features(ND = i)

    return t

## ETE3 TREE-VIZ FUNCTIONS ##

# basic tree style
tree_style = TreeStyle()
tree_style.show_leaf_name = True
tree_style.show_branch_length = False
tree_style.draw_guiding_lines = True
tree_style.complete_branch_lines_when_necessary = True

# make tree grow upward
#tree_style.rotation = 270
# and make it appear ultrametric (which it is!)
#tree_style.optimal_scale_level = "full"
tree_style.scale = 1200 # 3000 works nicely for full cteno tree
branchWidth = 5

# internal node style
nstyle = NodeStyle()
nstyle["fgcolor"] = "black"
nstyle["size"] = 0
nstyle["hz_line_width"] = branchWidth
nstyle["vt_line_width"] = branchWidth

# terminal node style
nstyle_L = NodeStyle()
nstyle_L["fgcolor"] = "black"
nstyle_L["size"] = 0
nstyle_L["hz_line_width"] = branchWidth
nstyle_L["vt_line_width"] = branchWidth

### Draw a tree with nodes color-coded and labeled by trait value
def traitTree(tree, traits, mapper, outDir, floatSwitch=0):
    ### Take dict of traits and [R,G,B]-returning function
    ### Draw a tree with the continuous trait painted on via a colormapping function

    def rgb2hex(r, g, b):
        hex = "#{:02x}{:02x}{:02x}".format(r, g, b)
        return hex

    for n in tree.traverse():
        if n.is_leaf():
            n.set_style(nstyle_L)
            # a manual way to set leaf labels
            #n.add_face(TextFace(str(n.name)), column=0, position="aligned")
        else:
            n.set_style(nstyle)
        #nd = TextFace(str(n.ND)) # label with node ID
        # make the face
        if round == 0:
            nd = TextFace(str(int(round(traits[n.ND])))) # label with rounded continuous trait value
        else:
            nd = CircleFace(10, color = rgb2hex(*[int(val) for val in mapper(traits[n.ND])]))
            #nd = TextFace(str(round(traits[n.ND], floatSwitch)))  # label with rounded continuous trait value

        #nd.background.color = rgb2hex(*[int(val) for val in mapper(traits[n.ND], gamma=0.8, scaleMax=255)]) # setup for wl2RGB
        #nd.background.color = rgb2hex(*[int(val) for val in mapper(traits[n.ND])])  # setup for grayscale
        nd.background.color = None
        #nd.margin_right = 2
        #nd.margin_top = 1
        #nd.margin_left = 2
        #nd.margin_bottom = 1
        nd.border.width = None
        nd.hz_align = 2
        #nd.inner_border.width = 1
        # outline for the circle node
        #ol = CircleFace(11, color = "black")
        #ol.hz_align = 0
        #n.add_face(ol, column=0, position="float-behind") #float-behind
        n.add_face(nd, column=0, position="float") #float
        # this is necessary for some reason to keep the tree from collapsing
        n.add_face(TextFace("       "), column=0, position="branch-bottom")

    outFile = outDir + "/cont_trait.pdf"
    tree.render(outFile, tree_style=tree_style)
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

### Draw test trees with the trait mapped to them and blocks on transition branches
def traitTestTrees(scenarios, traits, mapper, outDir, treePath, floatSwitch=0):

    # styles for convergent and transition nodes
    styleConv = NodeStyle()
    styleConv["fgcolor"] = "black"
    styleConv["size"] = 0
    styleConv["hz_line_width"] = branchWidth
    styleConv["vt_line_width"] = branchWidth
    styleConv["hz_line_color"] = "orange"
    styleConv["vt_line_color"] = "orange"

    styleTran = NodeStyle()
    styleTran["fgcolor"] = "black"
    styleTran["size"] = 0
    styleTran["hz_line_width"] = branchWidth
    styleTran["vt_line_width"] = branchWidth
    styleTran["vt_line_color"] = "orange"

    # face object indicating a transition
    block = RectFace(4, 20, "blue", "blue")
    block.margin_right = 2

    def rgb2hex(r, g, b):
        hex = "#{:02x}{:02x}{:02x}".format(r, g, b)
        return hex

    # loop thru trees
    for cutoff in sorted(scenarios.keys()):
        tree = init_tree(treePath)

        #tree.convert_to_ultrametric(tree_length=1)

        # parse scenario
        manual_mode_nodes = {"T": [], "C": []}
        p_events = scenarios[cutoff].strip().split("/")
        for e in p_events:
            l_e = map(int, e.split(","))
            manual_mode_nodes["T"].append(l_e[0])
            manual_mode_nodes["C"].extend(l_e[1:])

        # loop thru nodes
        for n in tree.traverse():

            #outline for the circle node
            #ol = CircleFace(11, color = "black")#
            #n.add_face(ol, column=1, position="float-behind") #float-behind

            # draw circle face with mapped color
            nd = CircleFace(14, color=rgb2hex(*[int(val) for val in mapper(traits[n.ND])]))
            nd.background.color = None
            nd.border.width = None

            if manual_mode_nodes:
                if n.ND in manual_mode_nodes["T"]:
                    n.set_style(styleTran)
                    n.add_face(block, column = 0, position="float")
                elif n.ND in manual_mode_nodes["C"]:
                    n.set_style(styleConv)
                else:
                    n.set_style(nstyle)
            else:
                nd.background.color = "white"

            n.add_face(nd, column=2, position="float")  # float
            # this is necessary for some reason to keep the tree from collapsing
            n.add_face(TextFace("       "), column=2, position="branch-bottom")
            n.img_style["size"] = 0

        outFile = str(round(cutoff, floatSwitch)).replace('.','_') + ".pdf"

        # prepend path to filename
        outFile = outDir + '/' + outFile
        tree.render(outFile, h=6, units="in", tree_style=tree_style)
        print >> sys.stderr, outFile
        # no return