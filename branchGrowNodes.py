# [nodes,neur_ves] = branchGrowNodes(nodes,neur_ves,params,idx,dir)
#
# Function to grow surface vasculature. Given a starting node (labeled as
# idx and dir), blood vessels are grown and branched with a distribution 
# given by the parameters set while avoiding crossover of vessels. The 
# inputs to this function are:
# 
# - nodes       - struct array containing blood vessel node information
#     .num      - identifying number
#     .root     - root node, if it exists (0 is no root)
#     .conn     - connecting node(s)
#     .pos      - position of node
#     .type     - string, type of node (edge, surf, vert, sfvt (surf/vert), capp) 
#     .misc     - used for misc. parameters
# - neur_ves    - surface vasculature occupance (for crossover avoidance)
# - params      - branching parameters
#     .maxit    - Maximum iteration to place nodes
#     .lensc    - Average distance between vasculature branch points
#     .varsc    - Standard deviation of distances between vascualure 
#                 branch points
#     .mindist  - Minimum inter-node distance
#     .varpos   - Standard deviation of vasculature placement
#     .dirvar   - The maximum branching angle
#     .branchp  - Probability of branching surface vasculature
#     .vesrad   - Radius of surface vasculature
# - idx         - starting branch node index
# - dir         - starting branch node growth direction
#
# And the outputs are
# - nodes       - Updated struct array containing blood vessel node information
#     .num      - identifying number
#     .root     - root node, if it exists (0 is no root)
#     .conn     - connecting node(s)
#     .pos      - position of node
#     .type     - string, type of node (edge, surf, vert, sfvt (surf/vert), capp) 
#     .misc     - used for misc. parameters
# - neur_ves    - Updated surface vasculature occupance (for crossover avoidance)
#
# 2017 - Alex Song

###########################################################################
import numpy as np

def branchGrowNodes(nodes,neur_ves,params,idx,dir):

    borderflag  = 1                                                           # Initialize flag to indicate being at a border
    overlapflag = 0                                                           # Choose whether or not to allow overlapping paths
    numIt       = 0                                                           #
    prevPos     = nodes(idx).pos(range(1,2))                                          #
    prevNum     = nodes(idx).num                                              #
    testIdxs2   = []                                                          #
    nv_size     = len(neur_ves)                                              #
    branchP     = 0                                                           #
    se          = strel('disk',params.vesrad)                                 #
    neur_ves2   = imdilate(neur_ves,se)                                       # Dilate the paths


    while(numIt<params.maxit && borderflag):
        if(np.random.rand<branchP):
            branchP = 0
            dirB    = dir - abs((0.5+np.random.rand)*params.dirvar)
            dir     = dir + abs((0.5+np.random.rand)*params.dirvar)
            [nodes,neur_ves] = branchGrowNodes(nodes,neur_ves,params,prevNum,dirB)
            overlapflag = 1
        else:
            branchP = branchP + params.branchp

        dirVect = [np.cos(dir), np.sin(dir)]
        vesDist = max(params.varsc*np.random.randn+params.lensc,params.mindist)
        nodePos = dirVect*vesDist+prevPos+params.varpos*np.random.randn(1,2)
        if(sum(nodePos<1) | sum(nodePos>nv_size(range(1,2)))):
            if(nodePos(1)<1):
                nodePos(1)=1
            if(nodePos(2)<1):
                nodePos(2)=1
            if(nodePos(1)>nv_size(1)):
                nodePos(1)=nv_size(1)
            if(nodePos(2)>nv_size(2)):
                nodePos(2)=nv_size(2)
            borderflag = 0
        testSubs = round(veclinspace(nodePos,prevPos,np.ceil(vesDist)))          # Sum of positions from linear interpolation to next position should not overlap
        testSubs = [testSubs bsxfun(@plus,testSubs,[01]) bsxfun(@plus, testSubs,[10])]
        testSubs(1,testSubs(1,:)>nv_size(1)) = nv_size(1)
        testSubs(2,testSubs(2,:)>nv_size(2)) = nv_size(2)
        testIdxs = sub2ind(nv_size(1:2),testSubs(1,:),testSubs(2,:))
        if(sum(neur_ves2(testIdxs))==0 || overlapflag):
            nodeNum        = len(nodes)+1
            nodes(nodeNum) = gennode(nodeNum,prevNum,prevNum,nodePos,'surf')
            prevPos        = nodePos
            prevNum        = nodeNum
            numIt          = numIt+1
            if(!isempty(testIdxs2)):
                neur_ves(testIdxs2) = 1
            testIdxs2 = testIdxs
        else:
            borderflag = 0
            testIdxs   = []
            testIdxs2  = []
        overlapflag = 0
    neur_ves(testIdxs2) = 1

    return [nodes,neur_ves]
###########################################################################
###########################################################################
