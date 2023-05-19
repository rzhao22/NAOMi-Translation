import numpy as np
import gennode as gennode
import branchGrowNodes as branchGrowNodes
import vessel_dijkstra as vessel_dijkstra
import pseudoRandSample2D as pseudoRandSample2D
import delnode as delnode

## imdilate, strel, and pos2dists functions in this file are not translated
def growMajorVessels(nv, node_p, vp):

# growMajorVessels(nv,np,vp)
# return nodes,nv
#
# Function to sample vasculature locations and setup nodes and connections of
# vasculature (including connections to penetrating vessels)
#
# - nv           - class instance for number of vessel parameters
#   .size        - size of volume [px]
#   .szum        - size of volume [um]
#   .nsource     - number of source vessels from volume surface
#   .ncapp       - number of capilliaries
#   .nvert       - number of penetrating vessels
#   .nnodes      - number of nodes
#   .nconn       - number of edges
#   .nvert_sum   - number of penetrating vessel connections to capilliaries
# - np           - node placement parameters
#   .maxit       - Maximum iteration to place nodes
#   .lensc       - average distance between vasculature branch points 
#   .varsc       - standard deviation of distances between vascualure 
#                  branch points  
#   .mindist     - minimum inter-node distance 
#   .varpos      - standard deviation of vasculature placement
#   .dirvar      - maximum branching angle
#   .branchp     - probability of branching surface vasculature
#   .vesrad      - radius of surface vasculature
# - vp           - class instance for vasculature parameters
#   .depth_surf  - Depth into tissue of surface vasculature
#   .sepweight   - Set weight that guides how far, on average, the
#                  vasculature nodes are placed
#   .mindists    - set minimum distance between nodes
#   .distWeightScale - scaling factor for weight of node distance scaling 
#   .randWeightScale - scaling factor for weight of nodes (variability)
#   .maxcappdist - maximum capilliary distance
#
# - nodes        - class instance containing blood vessel node information
#   .num         - identifying number
#   .root        - root node, if it exists (0 is no root)
#   .conn        - connecting node(s)
#   .pos         - position of node
#   .type        - string, type of node (edge, surf, vert, sfvt, capp) 
#   .misc        - used for misc. parameters
#
# 2017 - Alex Song

###########################################################################

    ## Initialize a few points for vertical vessels, Initialize some points in surface for surface vessels
    ## Create node class instance and grow initial (large) vessels from edges
    nodes = gennode()                               
    for i in range(nv.nsource):
        TMPIDX = (np.random.rand(1,2) >= [(nv.size[1]/(nv.size[0]+nv.size[1])), 0.5])
        if  (TMPIDX[0] and  TMPIDX[1]):
            TMPPOS = [np.ceil(nv.size[0]*np.ranodm.rand), 1, vp.depth_surf]
        elif  (TMPIDX[0] and not TMPIDX[1]):
            TMPPOS = [np.ceil(nv.size[0]*np.ranodm.rand), nv.size[1], vp.depth_surf]
        elif (not TMPIDX[0] and  TMPIDX[1]):
            TMPPOS = [1, np.ceil(nv.size[1]*np.ranodm.rand), vp.depth_surf]
        elif (not TMPIDX[0] and not TMPIDX[1]):
            TMPPOS = [nv.size[0],  np.ceil(nv.size[1]*np.ranodm.rand), vp.depth_surf]
        nodes[i] = gennode(i,0,[],TMPPOS,'edge',TMPIDX)

    neur_surf = False(nv.size[0:1])
    for i in range(nv.nsource):
        if (nodes[i].misc[0] and  nodes[i].misc[1]):
            randDir = 0.5*np.pi+np.random.randn*node_p.dirvar
        elif  (nodes[i].misc[0] and not nodes[i].misc[1]):
            randDir = 1.5*np.pi+np.ranodm.randn*node_p.dirvar
        elif (not nodes[i].misc[0] and  nodes[i].misc[1]):
            randDir = 0.0*np.pi+np.random.randn*node_p.dirvar
        elif (not nodes[i].misc[0] and not nodes[i].misc[1]):
            randDir = 1.0*np.pi+np.random.randn*node_p.dirvar
        [nodes,neur_surf] = branchGrowNodes(nodes,neur_surf,node_p,i,randDir)  
    nv.nlinks = len(nodes)-nv.nsource

    for i in range((nv.nlinks+nv.nsource)):
        i = i + nv.nsource
        nodes[i].pos = round([nodes[i].pos, vp.depth_surf])

    ## Sample additional locations such that diving vessels can semi-uniformly cover the surface
    neur_surf = imdilate(neur_surf,strel('disk',node_p.vesrad*2))

    surfpos = pseudoRandSample2D(nv.size[0:1],nv.nsurf,vp.mindists[0],vp.sepweight,(1-neur_surf).astype(np.float32)) # 
    surfpos = np.concatenate(surfpos,vp.depth_surf*np.ones((nv.nsurf,1)), axis=1)   
    surfpos = np.concatenate(np.reshape([nodes[0:nv.nlinks+nv.nsource].pos],[],nv.nlinks+nv.nsource).H,surfpos, axis=0)

    # Forms a block matrix of connections:
    # M = [A B]
    #     [C D]
    surfmat = pos2dists(surfpos)
    surfmat[0:nv.nlinks+nv.nsource,0:nv.nlinks+nv.nsource] = np.Inf
    surfmat[0:nv.nsource,0:nv.nsource] = 0
    for i in range(nv.nlinks+nv.nsource):
        if(nodes[i].root>0):
            surfmat[nodes[i].root,i] = np.norm(nodes[i].pos-nodes[nodes[i].root].pos)
            surfmat[i,nodes[i].root] = surfmat[nodes[i].root,i]
        
    surfmat[nv.nlinks+nv.nsource+1:,0:nv.nlinks+nv.nsource] = np.Inf

    TMPsurfmat = (surfmat**vp.distWeightScale)*(1+
                        vp.randWeightScale*np.random.randn(surfmat.shape))
    [_,surfpath] = vessel_dijkstra(TMPsurfmat,1)
    surfpath[0:nv.nsource]  = range(nv.nsource)
    for i in range(nv.nsurf+nv.nsource):
        if surfpath[i] == 0:
            [_,surfpath[i]] = min(surfmat[i,0:nv.nsource])
        
    

    for i in range(nv.nsurf):
        i = nv.nlinks+nv.nsource+i
        if (sum(surfpos[i,1]==[1, nv.size[0]]) or sum(surfpos[i,2]==[1, nv.size[1]])):
            nodes[i] = gennode(i,surfpath[i],surfpath[i],surfpos[i,:],'edge')    
        else:
            nodes[i] = gennode(i,surfpath[i],surfpath[i],surfpos[i,:],'surf')
    
    

    nv.nnodes = nv.nlinks+nv.nsource+nv.nsurf
    for i in range(nv.nnodes):
        if(nodes[i].root>0):
            nodes[nodes[i].root].conn = np.union(nodes[nodes[i].root].conn,i)
    

    ## Prune some surface vasculature and choose diving vessels
    se = strel('disk',round(vp.mindists[0]*2))
    neur_vert = False(nv.size[:1])
    for i in range(nv.nnodes):
        if ((nodes[i].type == 'surf') and (len(nodes[i].conn)==1)):
            if((neur_vert[nodes[i].pos[0],nodes[i].pos[1]]==0)):
                nodes[i].type = 'sfvt'
                TMP = False(nv.size[:1])
                TMP[nodes[i].pos[0],nodes[i].pos[1]] = 1
                neur_vert = neur_vert+imdilate(TMP,se)
            else:
                nodes = delnode(nodes,i)


    nodes = nodes[:nv.nnodes]
    surfidx = [i for i, node in enumerate(nodes) if node.type == 'sfvt']
    surfpos = np.reshape([nodes[(({nodes.type}=='surf'))].pos],3,[]).H
    surfpos = np.ravel_multi_index((surfpos[:,0],surfpos[:,1]), nv.size[:1])
    while (sum(node.type == 'sfvt' for node in nodes)< nv.nvert and sum(neur_vert(surfpos)==0)>0):
        TMPIDX = np.where(neur_vert(surfpos)==0)
        TMPIDX = surfidx(TMPIDX(np.ceil(np.random.rand*len(TMPIDX))))
        nodes(TMPIDX).type = 'sfvt'
        TMP = False(nv.size[:1])
        TMP[nodes[TMPIDX].pos[0],nodes[TMPIDX].pos[1]] = 1
        neur_vert = neur_vert+imdilate(TMP,se)

    ## Grow diving vessels to bottom of volume
    vertidx = [i for i, node in enumerate(nodes) if node.type == 'sfvt']
    TMPIDX = nv.nnodes
    for i in range(len(vertidx)):
        curr_node = vertidx[i]
        while(nodes[curr_node].pos[2]<nv.size[2]):
            TMPIDX = TMPIDX + 1
            node_pos = nodes[curr_node].pos+np.ceil([np.random.randn(1,2), 1]*[node_p.varpos, node_p.varpos, 
                     max(node_p.varsc*np.random.randn+node_p.lensc,node_p.mindist)])
            node_pos = min(max(node_pos,[1, 1, 1]),nv.size)
            nodes[TMPIDX] = gennode(TMPIDX,curr_node,curr_node,node_pos,'vert')
            nodes[nodes[TMPIDX].root].conn = np.union(nodes[nodes[TMPIDX].root].conn,TMPIDX)
            curr_node = TMPIDX

    nv.nvert = sum(node.type == 'sfvt' for node in nodes)
    nv.nvertconn = curr_node-nv.nnodes
    nv.nnodes = curr_node

    ## End nodes are initialized to a size
    ends = [i for i, node in enumerate(nodes) if len(node.conn) == 1]
    for i in range (len(ends)):
        nodes[ends[i]].misc = vp.vesSize[2]+np.random.gamma(3,(vp.vesSize[1]-vp.vesSize[2])/3) # gamma distribution with shape parameter 3 to set the vessel distribution

    return nodes, nv
