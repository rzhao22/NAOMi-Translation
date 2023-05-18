import scipy.sparse as sp

def nodesToConn(nodes):

# conn = nodesToConn(nodes)
#
# Function that generates the connection structure from the node structure
# for blood vessels. The input is:
#
# - nodes       - struct array containing blood vessel node information
#     .num      - identifying number
#     .root     - root node, if it exists (0 is no root)
#     .conn     - connecting node(s)
#     .pos      - position of node
#     .type     - string, type of node (edge, surf, vert, sfvt (surf/vert), capp) 
#     .misc     - used for misc. parameters
#
# The output is:
# - conn        - struct array containing blood vessel edge information
#     .start    - start node
#     .ends     - end node
#     .weight   - node weight
#     .locs     - position of locations
#     .misc     - used for misc. parameters
#
# 2017 - Alex Song

###########################################################################

    ends = [i for i, node in enumerate(nodes) if len(node.conn) == 1]
    connmat = sp.csr_matrix((length(nodes), length(nodes), dtype=float)

    for i = in range(len(ends)):
        curr_node = ends[i];
        while(nodes[curr_node].root>0):
            connmat[nodes[curr_node].num,nodes[curr_node].root) = ...
            np.sqrt(connmat[nodes[curr_node].num,nodes[curr_node].root)**2+nodes[ends[i]].misc**2)
            curr_node = nodes[curr_node].root

    # find the row and column indices of non-zero elements
    x, y = connmat.nonzero()
    # get non-zero values
    w = connmat[x, y].toarray().flatten

    conn = genconn()
    for i = in range(len(x)):
        conn[i] = genconn(x[i],y[i],w[i])
        
    return conn
