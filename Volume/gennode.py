import numpy as np

def gennode(num = [], root = [], conn = [], pos = [], type = [], misc = []):

# gennode(num,root,conn,pos,type,misc)
# return node
#
# Function to create a node class instance (branch point of blood vessels).  Nodes 
# are arranged as a class instance with fields: 
#
#   node.num      - identifying number
#   node.root     - root node, if it exists (0 is no root)
#   node.conn     - connecting node(s)
#   node.pos      - position of node
#   node.type     - string, type of node (edge, surf, vert, sfvt (surf/vert), capp) 
#   node.misc     - used for misc. parameters
#
# 2017 - Alex Song
#
###########################################################################
###########################################################################
## Set node values
    class Node:
        def __init__(self,num,root,conn,pos,type,misc):
            self.num  = num
            self.root = root
            self.conn = conn
            self.pos  = pos
            self.type = type
            self.misc = misc
    node = Node(num,root,conn,pos,type,misc)

    return node



###########################################################################
###########################################################################