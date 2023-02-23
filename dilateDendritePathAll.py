
def dilateDendritePathAll(paths,pathnums,obstruction):

# [paths,pathnums] = dilateDendritePathAll(paths,pathnums,obstruction)
#
# Function to simultaneously dilate all the dendrite paths in the volume.
# This function looks for empty space around the established dendrite
# positions to allocate to the expanded dendrite size, within taxicab 3
# The inputs of this function are:
#
# - paths       - Full set of simulated paths
# - pathnums    - Corresponding cell number
# - obstruction - Occupied space in volume
#
# The ouputs of this function are
#
# - paths       - Updated set of simulated paths
# - pathnums    - Updated corresponding cell number
#
# 2017 - Alex Song

###########################################################################

    # maximum number of voxels that an object can be dilated
    maxDist = 20

    meshrange = np.arange(-maxDist,maxDist+1)
    [x,y,z] = meshgrid(meshrange,meshrange,meshrange)
    dists = (np.power(x,2)+np.power(y,2)+np.power(z,2))
    dsz = dists.shape
    didx = np.argsort(dists.flatten())
    dval = np.sort(dists.flatten())
    dpos = np.nonzero(np.diff(dval))

    
    # python doesn't support single-precision       paths = single(paths);
    paths[obstruction!=0] = None
    dims = paths.shape
    pdims = np.prod(dims,axis=0)
    dshifts = np.array([-dims[0]*dims[1],dims[0]*dims[1],-dims[0],dims[0],-1,1]).T

    idxs = np.nonzero(paths>1)
    i = 0
    while(i<np.power(maxDist,2) and np.any(idxs))
        # shifts from center representing 
        [dx,dy,dz] = ind2sub(dsz,didx(dpos[i]+1:dpos[i+1])); #this needs to be fixed
        dx = dx-maxDist-1
        dy = dy-maxDist-1
        dz = dz-maxDist-1
        # shifts for tested indexesshifts for tested indexes
        jidxs = (dz)*dims[0]*dims[1]-(dy)*dims[0]-dx
        for j in range(0,len(idxs)):
            # make sure indexes are not out of bounds
            pidxs = idxs[j]+jidxs
            pidxs = pidxs[pidxs>0]
            pidxs = pidxs[pidxs<pdims]
            pidxs = pidxs[paths[pidxs]==0]
            didxt = np.zeros(len(pidxs)).T
            numval = pathnums[idxs[j]]
        # make sure that indexes are connected to rest of expanded process
        for k in range(len(pidxs)):
          didxs = pidxs[k]+dshifts
          didxs = didxs[didxs>0]
          didxs = didxs[didxs<pdims]
          if(sum(pathnums[didxs]==numval)):
            didxt[k] = True
        pidxs = pidxs[didxt]
        if(pidxs.any()):
          while(paths[idxs[j]]>1 and pidxs.any()):
            ridx = randi[len(pidxs)]
            pidx = pidxs[ridx]
            pidxs = np.delete(pidxs,ridx)
            paths[idxs[j]] = paths[idxs[j]]-1
            paths[pidx] = 1
            pathnums[pidx] = numval

      idxs = np.where(paths>1)
      i = i+1
      
    return paths, pathnums
