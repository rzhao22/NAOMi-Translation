import numpy as np

def multiBoundaryImage(X,varargin):

    X = np.array(X)
    
    if nargin > 1
        t = varargin{1}
    else
        t = 0.01


    I = cell(numel(X,3),1)

    for kk = 1:size(X,3)
        TMP       = np.array(X(:,:,kk))
        TMP       = TMP > t*max(TMP(:))
        if sum(TMP(:)) > 0
            [row,col] = find(TMP == 1,1,'first')
            I{kk}     = bwtraceboundary(TMP,[row, col],'N')
        else
            I{kk} = []
        
    
    return I

