import numpy as np
from numpy import fft
# masked_3DGP_test(grid_sz, l_scale, p_scale, mu, bin_mask, threshold, l_weights,type)
# return gp_tmp

def masked_3DGP_test(grid_sz, l_scale, p_scale, mu, bin_mask=1, threshold=1e-10, l_weights=1,type='single'):


    # [gp_vals] = masked_3DGP_v2(grid_sz, l_scale, p_scale, varargin)
    # 
    # This function draws from a 3-D GP with mean mu, length-scale l_scale and
    # variance scale parameter p_scale. Basically it returns a volume X, where
    #             X ~ N(mu, C), where C_{i,j} = p*e^{-(i-j)^2/(2*l^2)}
    # The method used is to draw i.i.d. random variables in the frequency
    # domain, apply the fft of C_{i,j} point-wise, and ise in ifft to return to
    # the spatial domain. This function was modified for speedups and multiple
    # spatial scales
    # 
    # It takes in parameters 
    #     - grid_sz:  The dimensions of the sample to take (if grid_sz is a
    #                 scalar, then the samples in each dimension will be that
    #                 value).
    #     - l_scale:  The length scale of the GP (i.e. l in the above
    #                 definition of C_{i,j})
    #     - p_scale:  The covariance scaling of the GP (i.e. p in the above
    #                 definition of C_{i,j})
    #     - mu:       The mean of the grid size. Can be either a scalar or a
    #                 matrix of the same size as X.
    #     - bin_mask: OPTIONAL Binary mask to set certain values of the output
    #                 to zero
    # 
    # And returns the output
    #     - X:        The sample from the GP, as defined above.
    #
    # 2017 - Adam Charles and Alex Song
    #
    ###########################################################################
    ## Check Inputs

    if not type:
        type = 'single'
    
    if not l_weights:
        l_weights = 1
    
    if not threshold:
        threshold = 1e-10
    
    if not bin_mask:
        bin_mask = 1
    
    if np.size(grid_sz) == 1:
        grid_sz = grid_sz*[1,1,1]

    if np.size(l_scale,1) == 1:
        l_scale = l_scale*[1,1,1]
    
    if np.size(l_weights) == 1:
        l_weights = np.ones([np.size(l_scale,0),1])*l_weights


    ###########################################################################
    ## Create Kernel

    wmx = np.pi/2
    if(type == 'single'):
        grid_x = np.reshape(np.square((np.linspace(-wmx,wmx,grid_sz[0])).astype('float32')),[-1,1,1], order='F')
        grid_y = np.reshape(np.square((np.linspace(-wmx,wmx,grid_sz[1])).astype('float32')),[1,-1,1], order='F')
        grid_z = np.reshape(np.square((np.linspace(-wmx,wmx,grid_sz[2])).astype('float32')),[1,1,-1], order='F')
    elif(type == 'double'):
        grid_x = np.reshape(np.square(np.linspace(-wmx,wmx,grid_sz[0])),[-1,1,1], order='F')
        grid_y = np.reshape(np.square(np.linspace(-wmx,wmx,grid_sz[1])),[1,-1,1], order='F')
        grid_z = np.reshape(np.square(np.linspace(-wmx,wmx,grid_sz[2])),[1,1,-1], order='F')
    
    
    # gp_vals = zeros(grid_sz,'single');
    gp_tmp = np.zeros(grid_sz,type)

    for i in range(np.size(l_scale,0)):
        ker_x = np.exp(np.square(-grid_x*l_scale(i,1)))
        ker_y = np.exp(np.square(-grid_y*l_scale(i,2)))
        ker_z = np.exp(np.square(-grid_z*l_scale(i,3)))
        gp_tmp = gp_tmp+np.prod(l_scale[i,:])*(np.square((ker_x * ker_y) * ker_z))*np.square(l_weights[i])
    
    gp_tmp = np.sqrt(gp_tmp)

    # if nargout > 1
    #     varargout{1} = gp_tmp;
    # end

    # sprev = rng
    # rng('shuffle','simdTwister')
    np.random.seed(10)
    for i in range(np.size(gp_tmp,2)):
        gp_tmp[:,:,i] = (np.random.randn(grid_sz[:2],type)+1j*np.random.randn(grid_sz[:2],type))*gp_tmp[:,:,i]
    
    # gp_tmp = exp(1j*2*pi*rand(grid_sz,type)).*gp_tmp;
    # gp_tmp = single(sqrt(chi2rnd(2,grid_sz)).*gp_tmp);
    # rng(sprev)
    np.random.seed(10)

    gp_tmp = 2*p_scale*bin_mask*np.sqrt(np.size(gp_tmp))*np.real(fft.ifftshift(fft.ifftn(fft.ifftshift(gp_tmp))))+mu # Move random process back to the spatial domain

    return gp_tmp

    ###########################################################################
    ###########################################################################

    #     ker_1 = bsxfun(@times,bsxfun(@times,ker_x,ker_y),ker_z);
    #     all_ker = all_ker+prod(l_scale(i,:))*ker_1.^2*l_weights(i).^2;

    # gp_vals = randn(grid_sz,type)+1j*randn(grid_sz,type);              # First make a random array
    # gp_vals = gp_vals.*all_ker;

    # gp_vals = p_scale*(2^4.5/pi^1.5)*bin_mask.*gp_vals/sqrt(length(l_weights))+mu;  # Apply a binary mask as needed