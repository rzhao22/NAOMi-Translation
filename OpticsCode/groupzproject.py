import numpy as np
def groupzproject(images, groupsize, type = 'mean'):
# groupzproject(images,groupsize,type)
# return out
#
# Function to project a 3D matrix along the third dimension in groups of
# the specified size. The inputs are
#
#   - images      - 3D matrix to project in groups
#   - groupsize   - Size of projection group. Needs to be a factor of
#                   size(images,3)
#   - type        - One of 'sum' 'prod' 'mean' 'median' 'mode' 'max' 'min' 
#                   for projection type. Default is 'mean'
#
# The output is
#   - out         - group projected 3D matrix
#
# 2017 - Alex Song

###########################################################################

    # if nargin < 3 type = 'mean' end                                          # If the projection type is unspecified, use the mean

    if(np.size(images,2) % groupsize==0):
        if type == 'sum':
            out = np.squeeze(np.sum(np.reshape(images,(np.size(images,0),np.size(images,1),
                groupsize,np.size(images,2)/groupsize)),axis = 2))
        elif type == 'prod':
            out = np.squeeze(np.prod(np.reshape(images,(np.size(images,0),np.size(images,1), 
                groupsize,np.size(images,2)/groupsize)),axis=2))
        elif type == 'mean':
            out = np.squeeze(np.mean(np.reshape(images,(np.size(images,0),np.size(images,1), 
                groupsize,np.size(images,2)/groupsize)),axis=2))
        elif type == 'median':
            out = np.squeeze(np.median(np.reshape(images,(np.size(images,0),np.size(images,1), 
                groupsize,np.size(images,2)/groupsize)),axis = 2))
        elif type == 'mode':
            out = np.squeeze(np.mode(np.reshape(images,(np.size(images,0),np.size(images,1), 
                groupsize,np.size(images,2)/groupsize)),axis = 2))
        elif type == 'max':
            out = np.squeeze(np.amax(np.reshape(images,(np.size(images,0),np.size(images,1), 
                groupsize,np.size(images,2)/groupsize)), axis = 2))
        elif type == 'min':
            out = np.squeeze(np.amin(np.reshape(images,(np.size(images,0),np.size(images,1), 
                groupsize,np.size(images,2)/groupsize)), axis = 2))
        else:
            TypeError
            print('need a valid type')

    else:
        if type == 'sum':
            out = np.squeeze(np.sum(np.reshape(images[:,:,:groupsize*np.floor(np.size(images,2)/groupsize)], 
                (np.size(images,0),np.size(images,1), groupsize,np.floor(np.size(images,2)/groupsize))),axis = 2))
            out = np.concatenate((out,np.sum(images[:,:,np.size(images,3)-(np.size(images,2)%groupsize)+1:],axis = 2)), axis = 2)
        elif type == 'prod':
            out = np.squeeze(np.prod(np.reshape(images[:,:,:groupsize*np.floor(np.size(images,2)/groupsize)], 
                (np.size(images,0),np.size(images,1), groupsize,np.floor(np.size(images,2)/groupsize))),axis = 2))
            out = np.concatenate((out,np.prod(images[:,:,np.size(images,2)-(np.size(images,2)%groupsize)+1:],axis=2)), axis=2)
        elif type == 'mean':
            out = np.squeeze(np.mean(np.reshape(images[:,:,:groupsize*np.floor(np.size(images,2)/groupsize)], 
                (np.size(images,0),np.size(images,1), groupsize,np.floor(np.size(images,2)/groupsize))),axis=2))
            out = np.concatenate((out,np.mean(images[:,:,np.size(images,2)-(np.size(images,2)%groupsize)+1:],axis=2)), axis=2)
        elif type == 'median':
            out = np.squeeze(np.median(np.reshape(images[:,:,:groupsize*np.floor(np.size(images,2)/groupsize)], 
                (np.size(images,0),np.size(images,1), groupsize,np.floor(np.size(images,2)/groupsize))),axis=2))
            out = np.concatenate((out,np.median(images[:,:,np.size(images,2)-(np.size(images,2)%groupsize)+1:],axis=2)), axis=2)
        elif type == 'mode':
            out = np.squeeze(np.mode(np.reshape(images[:,:,:groupsize*np.floor(np.size(images,2)/groupsize)], 
                (np.size(images,0),np.size(images,1), groupsize,np.floor(np.size(images,2)/groupsize))),axis=2))
            out = np.concatenate((out,np.mode(images[:,:,np.size(images,2)-(np.size(images,2)%groupsize)+1:],axis=2)), axis=2)
        elif type == 'max':
            out = np.squeeze(np.amax(np.reshape(images[:,:,:groupsize*np.floor(np.size(images,2)/groupsize)], 
                (np.size(images,0),np.size(images,1), groupsize,np.floor(np.size(images,2)/groupsize))),axis=2))
            out = np.concatenate((out,np.amax(images[:,:,np.size(images,2)-(np.size(images,2)%groupsize)+1:],axis=2)), axis=2)
        elif type == 'min':
            out = np.squeeze(np.amin(np.reshape(images[:,:,:groupsize*np.floor(np.size(images,2)/groupsize)], 
                (np.size(images,0),np.size(images,1), groupsize,np.floor(np.size(images,2)/groupsize))),axis=2))
            out = np.concatenate((out,np.amin(images[:,:,np.size(images,2)-(np.size(images,2)%groupsize)+1:],axis=2)), axis=2)
        else:
            TypeError('need a valid type')
    return out

    ###########################################################################
    ###########################################################################