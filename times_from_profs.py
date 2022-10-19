def times_from_profs(mov, NeurProf, lambda, varargin):  

    if isempty(lambda):
        lambda = 0                                                            # Default to non-negative least-squares

    if nargin > 3:
        BGprofs = varargin{1}                                                 # Set background 
    else:
        BGprofs = []                                                          # Set background ROIs to null

    if nargin > 4:
        batch_size = varargin{2}                                              # Set batch size
        if isempty(batch_size):
            batch_size = 3000                                                 # default batch size of 3000
    else:
        batch_size = 3000                                                     # default batch size of 3000

    if nargin > 5:
        tfocs_Tol = varargin{3}                                               # Set batch size
        if isempty(tfocs_Tol):
            tfocs_Tol = 1e-10*batch_size*size(NeurProf,3)                     # default tolerance based on number of variables
    else:
        tfocs_Tol = 1e-10*batch_size*size(NeurProf,3)                         # default tolerance based on number of variables

    if nargin > 6:
        bin_mask = varargin{4}                                                # Allow for an optional binary mask (e.g. if some elements are not to be esimated)
    else:
        bin_mask = 1                                                          # Default is no mask


    if isempty(BGprofs):
        BGprofs = median(mov,3)                                               # If no backgroun ROIS, calculate median
        BGprofs = BGprofs/norm(BGprofs(:))
    else if isnan(BGprofs):
        BGprofs = []                                                          # If no backgroun ROIS, calculate median

    n_t   = size(mov,3)                                                       # Get number of frames in mov
    n_BGs = size(BGprofs,3)                                                   # Get number of background frames
    NeurProf_TMP = [reshape(NeurProf, [], size(NeurProf, 3)),
                                             reshape(BGprofs,[],n_BGs)]       # Reshape current pairs into matrix

    x0 = max(NeurProf_TMP.'*reshape(mov, [], n_t),0)                          # Set up initial projection

    ls_opt.tol        = tfocs_Tol                                             # Set tolerance parameter for TFOCS (1e-3 is standard)
    ls_opt.nonneg     = true                                                  # Set solution to be non-negative
    ls_opt.printEvery = 0                                                     # No outputs!

    ###########################################################################
    ## Run Optimization

    N_gaps = ceil(size(mov,3)/batch_size)                                     # Calculate number of batches
    T_gaps = round(linspace(0,size(mov,3),N_gaps+1))                          # Get batch limits

    if max(lambda(:)) < 0:
        # x_est = NeurProf_TMP\reshape(mov, [], n_t)      
        for kk in range(1, N_gaps)
            print('Performing Batch %d of %d.\n', kk, N_gaps)
            n_k = T_gaps(kk+1)-T_gaps(kk)                                     # Calculate the number of frames in this batch
            x_est(:,(T_gaps(kk)+1):T_gaps(kk+1)) = linsolve(NeurProf_TMP,
                      reshape(mov(:,:,(T_gaps(kk)+1):T_gaps(kk+1)), [], n_k)) # Use MATLAB's built-in solver to solve the non-negative least-squares problem for each movie frame.      

    else if max(lambda(:)) == 0
        for kk in range(1, N_gaps)
            print('Performing Batch %d of %d.\n', kk, N_gaps)
            n_k = T_gaps(kk+1)-T_gaps(kk)                                     # Calculate the number of frames in this batch
            if numel(bin_mask) == numel(mov):
                A_ls1 = linop_matrix(NeurProf_TMP, 'R2R', n_k)                # Set up linear operator as a matrix multiplication
                A_idx = linop_subsample(reshape(bin_mask(:,:,(T_gaps(kk)+1):
                                                      T_gaps(kk+1)), [], n_k))# Set up subsample
                A_ls  = linop_compose(A_idx,A_ls1)                            # Compose subsampling with matrix multiplication
            else:
                A_ls = linop_matrix(NeurProf_TMP, 'R2R', n_k)                 # Set up linear operator as a matrix multiplication
            x_est(:,(T_gaps(kk)+1):T_gaps(kk+1)) = solver_NNLS(A_ls, 
                reshape(mov(:,:,(T_gaps(kk)+1):T_gaps(kk+1)), [], n_k), 
                x0(:,(T_gaps(kk)+1):T_gaps(kk+1)), ls_opt)                    # Run NNLS via TFOCS: Call to TFOCS
        end
    else
        for kk in range(1, N_gaps)
            print('Performing Batch %d of %d.\n', kk, N_gaps)
            n_k = T_gaps(kk+1)-T_gaps(kk)                                     # Calculate the number of frames in this batch
            if numel(bin_mask) == numel(mov):
                A_ls1 = linop_matrix(NeurProf_TMP, 'R2R', n_k)                # Set up linear operator as a matrix multiplication
                A_idx = linop_subsample(reshape(bin_mask(:,:,(T_gaps(kk)+1):
                                           T_gaps(kk+1)), [], n_k))           # Set up subsample
                A_ls  = linop_compose(A_idx,A_ls1)                            # Compose subsampling with matrix multiplication
            else:
                A_ls = linop_matrix(NeurProf_TMP, 'R2R', n_k)                 # Set up linear operator as a matrix multiplication

            if numel(lambda)==1:
                lambda2 = [lambda*ones(size(NeurProf,3),1) zeros(n_BGs,1)]   # Give the median a zero value
                lambda2 = single(single(lambda2)*single(ones(1,n_k)))         # Make lambda the same size as x_all expects to be
            else if numel(lambda) == size(NeurProf, 3):                          # Case where each profile has it's own lambda
                lambda2 = [lambda(:) zeros(n_BGs,1)]                         # Give the median a zero value
                single(single(lambda2)*single(ones(1,n_k)))                   # Make lambda the same size as x_all expects to be
            else if numel(lambda) == n_t*size(NeurProf, 3):                      # Case where each profile and each time has it's own lambda
                lambda2 = [single(reshape(lambda(:,(T_gaps(kk)+1):T_gaps(kk+1)),
                    [],n_k))zeros(n_BGs,n_k,'single')]                       # Make lambda the same size as x_all expects to be
            else if numel(lambda) == (size(NeurProf, 3)+n_BGs):                  # Case where each profile and each background componant has it's own lambda
                lambda2 = single(single(lambda(:))*single(ones(1,n_k)))       # Make lambda the same size as x_all expects to be
            else if numel(lambda) == n_t*(size(NeurProf, 3)+n_BGs):              # Case where each profile and each background componant has it's own lambda at every time
                lambda2 = single(reshape(lambda(:,(T_gaps(kk)+1):T_gaps(kk+1)),
                                               [],n_k))                       # Make lambda the same size as x_all expects to be

            x_est(:,(T_gaps(kk)+1):T_gaps(kk+1)) = solver_L1RLS(A_ls, 
                reshape(mov(:,:,(T_gaps(kk)+1):T_gaps(kk+1)), [], n_k), lambda2, 
                x0(:,(T_gaps(kk)+1):T_gaps(kk+1)), ls_opt)                    # Run Lasso via TFOCS: Call to TFOCS. lambda should be the same size as the initial point x0

    ###########################################################################
    ## Output Parsing

    if nargout == 0:
        pass
    else if nargout == 1:
        varargout{1} = x_est                                                  # If one output, lump all time traces together
    else:
        varargout{1} = x_est(1:(end-n_BGs),:)                                 # If one output, separate out median trace
        varargout{2} = x_est((end-n_BGs+1):end,:)                                           

    return varargout

