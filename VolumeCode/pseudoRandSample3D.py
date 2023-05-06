import numpy as np

def pseudoRandSample3D(sz, nsamps, pdf, width = 2, weight = 1):
## function [pos,pdf] = pseudoRandSample3D(sz, nsamps, width, weight, pdf)

# [pos,pdf] = pseudoRandSample3D(sz, nsamps, width, weight, pdf)
# 
# Function to sample pseudo randomly within a 3D array, with partial
# exclusion of previously sampled locations. The inputs are:
#
#  - sz         - Size of input matrix
#  - nsamps     - Number of locations to sample
#  - width      - Width of Gaussian exclusionary zone
#  - weight     - Weight of Gaussian (0-1, 1 fully excludes previously
#                 sampled locations)
#  - pdf        - Probability density function of allowable locations to
#                 sample
# 
# The outputs are:
#  - pos        - Sampled positions
#  - pdf        - Resulting probability density function of future
#                 positions that may be sampled
# 
# 2017 - Alex Song

###########################################################################
    # if(nargin<5):
    #     pdf = np.ones(sz,'single')
    if not pdf:
        pdf = np.ones(sz).astype(np.float32)


    # if(nargin<4)
    # weight = 1
    # end

    # if(nargin<3)
    # width = 2
    # end

    X,Y,Z = np.meshgrid(np.linspace(-np.ceil(2*width),np.ceil(2*width)),np.linspace(-np.ceil(2*width),np.ceil(2*width)))
    gpdf = (-weight*np.exp(-(X**2+Y**2+Z**2)/(width**2))).astype(np.float32)
    pos = np.zeros((nsamps,3))
    i = 0
    while(i<nsamps):
        rndpt = np.ceil(np.random.rand(1,3)*sz)
        if((pdf(rndpt[0],rndpt[1],rndpt[2])-np.random.rand)>0):
            xc = np.int32([max(0,2*width+1-rndpt[0]), min(0,sz[0]-rndpt[0]-2*width)+4*width]+1)
            yc = np.int32([max(0,2*width+1-rndpt[1]), min(0,sz[1]-rndpt[1]-2*width)+4*width]+1)
            zc = np.int32([max(0,2*width+1-rndpt[2]), min(0,sz[2]-rndpt[2]-2*width)+4*width]+1)
            
            xi = np.int32([max(1,rndpt[0]-2*width), min(sz[0],rndpt[0]+2*width)])
            yi = np.int32([max(1,rndpt[1]-2*width), min(sz[1],rndpt[1]+2*width)])
            zi = np.int32([max(1,rndpt[2]-2*width), min(sz[2],rndpt[2]+2*width)])
            pdf[xi[0]:xi[1],yi[0]:yi[1],zi[0]:zi[1]] = pdf[xi[0]:xi[1],yi[0]:yi[1],zi[0]:zi[1]]+gpdf[xc[0]:xc[1],yc[0]:yc[1],zc[0]:zc[1]]
            pos[i,:] = rndpt
            i = i+1