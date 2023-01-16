# conn = genconn(start,ends,weight,locs,misc)
#
# Function to generate a connection structure (piece of blood vessel 
# edge between nodes) Connections are a structure array with fields: 
#
#   conn.start      - start node
#   conn.ends       - end node
#   conn.weight     - node weight
#   conn.locs       - position of locations
#   conn.misc       - used for misc. parameters
#
# 2017 - Alex Song

###########################################################################

def genconn(start = None, ends = None, weight = None, locs = None, misc = None):

    conn = {}

    nargin = 0

    if(start != None):
        nargin += 1
    if (ends != None):
        nargin += 1
    if (weight != None):
        nargin += 1
    if (locs != None):
        nargin += 1
    if (misc != None):
        nargin += 1

    if(nargin<1):
        start  = []
    if(nargin<2):
        ends   = []
    if(nargin<3):
        weight = []
    if(nargin<4):
        locs   = []
    if(nargin<5):
        misc   = []
    
    conn[start]  = start
    conn[ends]   = ends
    conn[weight] = weight
    conn[locs]   = locs
    conn[misc]   = misc
    return conn
