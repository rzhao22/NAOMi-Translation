import numpy as np
from dendrite_dijkstra_cpp import dendrite_dijkstra_cpp
##function [distance,pathfrom] = dendrite_dijkstra2(M,dims,root)
def dendrite_dijkstra2(M, dims, root):

# function [distance,pathfrom2] = dendrite_dijkstra2(M,dims,root)
# 
# Function to run Dijkstra's algorithm for growing out dendrites. Inputs to
# the function are
#    - M    - dims-size matrix that indicates blockages in the volume (i.e.
#             other dendrites, blood vessels, etc.)
#    - dims - 3x1 Size of the volume to grow dendrites in
#    - root - 3x1 Starting point of the path
#
# Outputs are:
#    - distance - Distance the path has traveled (dendrite length)
#    - pathfrom - Path through the volume that the dendrite takes
# 
# 2017 - Alex Song
#
###########################################################################

    try:
        root2 = np.ravel_multi_index((root[0], root[1], root[2]), dims)  #sub2ind(dims,root(1),root(2),root(3))                         # Get index set location of the root node from the subscript location provided
    except:
        print(root)
        raise ValueError('root')

    if(type(M) != type(np.float32)):
        M = (M).astype(np.float32)                              # Make sure M is a single
        raise TypeError('M must be a single, casting as a single')
        

    e  = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]])                        # Set the adjacent edges to loop through: R L U D T B (right/left/up/down/towards/back)
    pe = e*np.array([[1],[dims[0]],[dims[0]*dims[1]]])                                        # Set distance based on dimensions of the volumes

    [distance,pathfrom_tmp] = dendrite_dijkstra_cpp(M,pe.astype(np.int32),root2.astype(np.int32)) # Run the mexed dijkstra algorithm for speed

    distance  = np.reshape(distance,dims)                                        # Reshape outputs to the size of the volume
    pathfrom  = np.zeros((np.prod(dims),3))                                           # 
    pidxs     = np.where(pathfrom_tmp>0)                                              # 
    [pathfrom[pidxs,0],pathfrom[pidxs,1],pathfrom[pidxs,2]] = np.unravel_index(pathfrom_tmp[pidxs], dims) # Change the path locations from an index set to subscript indexing
    pathfrom = np.reshape(pathfrom,[dims, 3])

    return distance, pathfrom

###########################################################################
###########################################################################
