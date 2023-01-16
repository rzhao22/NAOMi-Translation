# nodes = delnode(nodes,num)
#
# Function to delete a single node (branch point of blood vessels). Nodes 
# are arranged as a structure array with fields:
#
#   node.num      - identifying number
#   node.root     - root node, if it exists (0 is no root)
#   node.conn     - connecting node(s)
#   node.pos      - position of node
#   node.type     - string, type of node (edge, surf, vert, sfvt (surf/vert), capp)
#   node.misc     - used for misc. parameters
#
# 2017 - Alex Song

###########################################################################

def delnode(nodes,num):
    # Iterate through all the nodes
    for i in nodes(num).conn:
        if(i == nodes(num).conn):
            for j in nodes(i).conn(num):
                if(nodes(i).conn==num):
                    j = [];
            if(nodes(i).root == num):
                nodes(i).root = []

    node = gennode()
    nodes(num).num = node.num
    nodes(num).root = node.root
    nodes(num).conn = node.conn
    nodes(num).pos = node.pos
    nodes(num).type = node.type
    nodes(num).misc = node.misc

    return nodes

###########################################################################
###########################################################################