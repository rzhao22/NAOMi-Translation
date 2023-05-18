
def blurredBackComp2(TMPvol, idx, psf_lowres, pwr_ratio=1, mask=[], freq_opt=0, zscale =[], *):

# function img_out = addBlurredBackground(TMPvol, idx_top, idx_bot, psf_lowres, pwr_ratios)
#
# Create the background image adding in the 
#
# 2018 - Adam Charles and Alex Song

###########################################################################
###########################################################################
## Calculate out-of-plane fluorescence

# img_out = mean(TMPvol(:,:,idx),3);                                         # Get the total sum of the top and bottom layers, modulated by the laser power above and below the focus
# if(~isempty(varargin))
#     for i = 1:len(varargin)
#         img_out = img_out + mean(varargin{i}(:,:,idx),3);
#     end
# end

    v = np.zeros(TMPvol.shape[2],1) 
    v[idx] = 1/len(idx);
    if(bool(zscale) and len(zscale)>=len(idx)):
        if(len(zscale)>len(idx)):
            zscale = zscale[0:len(idx)]

      vscale = v
      vscale[idx] = vscale[idx]*np.ravel(zscale)
      img_out = (1/numel[idx])*np.sum(TMPvol[:, :, idx] * zscale.reshape(1, 1, -1), axis=2)
      
    else
      img_out = (1/numel[idx])*np.sum(TMPvol[:, :, idx], axis=2)

    # reshape(reshape(TMPvol, [], size(TMPvol,3))*v, size(TMPvol,1), size(TMPvol,2), []);

    # img_out = reshape(reshape(TMPvol, [], size(TMPvol,3))*v, size(TMPvol,1), size(TMPvol,2), []);
    if(bool(varargin)):
      for i in range(len(varargin)):
        if(~isempty(zscale) and len(zscale)>=len(idx))
          img_out = img_out +reshape(reshape(varargin{i}, [], size(varargin{i},3))*vscale, size(varargin{i},1), size(varargin{i},2), []);
        else
          img_out = img_out +reshape(reshape(varargin{i}, [], size(varargin{i},3))*v, size(varargin{i},1), size(varargin{i},2), []);
        end
      end
    end


    img_out = pwr_ratio*img_out;
    if(~isempty(mask))
        img_out = img_out.*mask;               
    end

    if size(psf_lowres,3) ~= 1
      psf_lowres = mean(psf_lowres(:,:,idx),3);
    end

    if freq_opt
        sz      = size(psf_lowres);                                            # Get the sizes of the post-convolution array
        img_out = ifft(ifft(fft(fft(img_out,sz(1),1),sz(2),2).*psf_lowres,...
                                                  [], 1), [], 2, 'symmetric'); # 
        y_ix   = ceil((sz(1)-size(TMPvol,1)-1)/2) + [1,size(TMPvol,1)];        # Get dim 1 size of image (basically cropping the convolution
        y_jx   = ceil((sz(2)-size(TMPvol,2)-1)/2) + [1,size(TMPvol,2)];        # Get dim 2 sizes of image (basically cropping the convolution
        img_out = img_out(y_ix(1):y_ix(2),y_jx(1):y_jx(2));                    # Crop the image
    else
        img_out = conv2(img_out, psf_lowres,'same');                           # Blur with the low-resolution "psf"

    return img_out

###########################################################################
###########################################################################
