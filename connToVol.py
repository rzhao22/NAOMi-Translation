import random.randint as randint
import cv2

def connToVol(nodes, conn, nv, idxs=np.arange(0,len(conn)), neur_ves=np.zeros(np.shape(nv))):

# [neur_ves,conn] = connToVol(nodes,conn,nv,idxs,neur_ves)
# 
# Function to add vasculature connectivity to a volume. Inputs to this function
# are:
# 
# - nodes       - struct array containing blood vessel node information
#     .num      - identifying number
#     .root     - root node, if it exists (0 is no root)
#     .conn     - connecting node(s)
#     .pos      - position of node
#     .type     - string, type of node (edge, surf, vert, sfvt (surf/vert), capp) 
#     .misc     - used for misc. parameters
# - conn        - struct array containing blood vessel edge information
#     .start    - start node
#     .ends     - end node
#     .weight   - node weight
#     .locs     - position of locations
#     .misc     - used for misc. parameters
# - nv.size     - volume size (pixels)
# - idxs        - indexes of vessels to generate
# - neur_ves    - simulated blood vessel volume
# 
# The ouputs of this function are
# - neur_ves    - Updated simulated blood vessel volume
# - conn        - Updated struct array containing blood vessel edge information
#     .start    - start node
#     .ends     - end node
#     .weight   - node weight
#     .locs     - position of locations
# 
# 2017 - Alex Song

###########################################################################


    for j in range(0,len(idxs)):
        i = idxs(j)
        TMPL = np.setxor1d(nodes[conn[i].start].conn,conn[i].ends)
        if (TMPL.any()):
            TMPL = TPML[randint(len(TMPL)] #this doesn't rly make sense
        TMPU = np.setxor1d(nodes[conn[i].ends].conn,conn[i].start)
        if (TMPU.any()):
            TMPU = TPMU[randint(len(TMPU)] #this doesn't rly make sense

        numsamp = 2*np.linalg.norm(nodes[conn[i].start]).pos.T-nodes(conn(i).ends).pos.T) #need help with .pos and fnval
        if (TMPU.any()) and (TMPL.any()):
            spts = cscvn([nodes[TMPL].pos.T, nodes(conn[i].start).pos.T, nodes[conn[i].ends]).pos.T, nodes(TMPU).pos.T])
            ves_loc = np.ceil(fnval(spts,np.linspace(spts.breaks(2),spts.breaks(3),numsamp)).T)
        elif (TMPL.any()):
            spts = cscvn([nodes[TMPL].pos.T, nodes(conn[i].start).pos.T, nodes[conn[i].ends).pos.T])
            ves_loc = np.ceil(fnval(spts,np.linspace(spts.breaks(2),spts.breaks(3),numsamp)).T)
        elif(TMPU.any()):
            spts = cscvn([nodes[conn[i].start).pos.T, nodes(conn[i].ends).pos.T, nodes[TMPU].pos.T])
            ves_loc = np.ceil(fnval(spts,np.linspace(spts.breaks(1),spts.breaks(2),numsamp)).T)
        else:
            spts = cscvn([nodes[conn[i].start).pos.T, nodes(conn[i].ends).pos.T])
            ves_loc = np.ceil(fnval(spts,np.linspace(spts.breaks(1),spts.breaks(2),numsamp)).T)

                           
        ves_loc = np.maximum(ves_loc,[1, 1, 1])
        ves_loc = np.minimum(ves_loc,nv.size)
        ves_loc = np.unique(ves_loc,axis=0)
        conn[i].locs = ves_loc


        min_idx = np.maximum(ves_loc.min(axis=0)-np.ceil(conn[i].weight),[1, 1, 1])
        max_idx = np.minimum(ves_loc.max(axis=0)+np.ceil(conn[i].weight),nv.size)
        ves_loc = ves_loc-min_idx
  
        TMP = np.zeros(((max_idx-min_idx+1)[2],(max_idx-min_idx+1)[0],(max_idx-min_idx+1)[1]))

        for each in range(len(ves_loc[:,2])):
            TMP[ves_loc[:,2][each]][ves_loc[:,0][each],ves_loc[:,1][each]]=1
        
        [x,y,z] = np.meshgrid(-np.ceil(conn[i].weight):ceil(conn[i].weight)) #need help with the rest of this
        se      = strel(sqrt(x.^2 + y.^2 + z.^2) <=conn(i).weight) 
        TMP     = imdilate(TMP,se)
        neur_ves(min_idx(1):max_idx(1),min_idx(2):max_idx(2),min_idx(3):max_idx(3)) = ...
        neur_ves(min_idx(1):max_idx(1),min_idx(2):max_idx(2),min_idx(3):max_idx(3))+TMP

              

    return neur_ves, conn
