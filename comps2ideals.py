import numpy as np

def comps2ideals(comps, baseim, varargin):
    ideal 
    if nargin > 2 :
        k = varargin{1}
    else 
        k = 2
    
    minNumEl = 5            # Set the minimum number of
    compsz = size(comps)    # Get size of components
    if(length(compsz)==2):
        compsz = [compsz 1]
        
    ideal = numpy.array(np.reshape(bsxfun(@rdivide,comps,baseim),(-1,compsz(3))))
    
    
    for i in range(1, size(comps,3)):
        sort_i = np.sort(ideal(:,i)).copy()
        sort_i = sort_i[::-1]
        sort_i = sort_i[not isnan(sort_i)];
        cutoff = 1/(k+2/np.mean(sort_i(1 : minNumEl)));

        if(sum(ideal(:,i)>cutoff)):
            rp = regionprops(np.reshape(ideal(:,i) > cutoff, (compsz(1), -1)),'PixelIdxList')
            if(not isempty(rp)):
                max_array = np.array(max(cellfun(@length,{rp.PixelIdxList})))
                val = max_array[:,1]
                temp = np.zeros(size(ideal,1),1)
                if(length(rp(val).PixelIdxList) >= minNumEl):
                    temp(rp(val).PixelIdxList) = ideal(rp(val).PixelIdxList,i)

            ideal(:,i) = temp

        else
            ideal(:,i) = zeros(size(ideal,1),1)


    ideal = reshape(ideal,compsz)   # Reshape ideal components back to the sizes of the components
    ideal = (ideal > 0) * comps

    return ideal 