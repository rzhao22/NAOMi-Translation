# [mask,psfs3,psfs] = genCorticalLightPath(vol_params,psf_params,vol_out,Uin)
#
# This function generates a point-by-point map of obscuration for TPM
# imaging setting across a neural volume. The function generates a 3D mask
# that modulates the intensity along the light paths. The inputs are
# 
#   - vol_params      - Struct with parameters for the volume generation
#       .vol_sz       - 3-element vector with the len (in um) of the 
#                       volume to generate (default = 100x100x30um)
#       .vres         - resolution to simulate volume at (default = 2
#                       samples/um)
#       .vol_depth    - Depth of the volume under the brain surface
#   - psf_params      - Struct contaning the parameters for the PSF
#       .n_diff       - Shift in refractive index from vessels to tissue
#       .lambda       - Two-photon excitation wavelength (um)
#       .obj_fl       - Objective focal length (mm)
#       .ss           - Subsampling factor for fresnel propagation
#       .sampling     - Spatial sampling for tissue occlusion mask
#       .psf_sz       - Default two-photon PSF len simulated (um)
#       .prop_sz      - Fresnel propagation length outside of volume (um)
#   - vol             - Simulated volume impacting propagation (vessels)
#   - Uin             - Input scalar field
#
# The outputs are
#   - mask            - 2D mask giving relative two-photon excitation at
#                       each position laterally
#   - psfs3           - Average aberrated PSF across the simulated field
#   - psfs            - All PSFs at each simulated position
#
# 2016 - Alex Song and Adam Charles

###########################################################################
## Parameters for light path obscurations
import numpy as np

def genCorticalLightPathLite(vol_params,psf_params,phzA,phzB,phzC,Uin):

  '''
  if(nargout>5):
    psf2p2 = [];
  '''

  vres      = vol_params.vres;
  vol_sz    = vol_params.vol_sz;
  vol_depth = vol_params.vol_depth;
  vasc_sz   = vol_params.vasc_sz;
  verbose   = vol_params.verbose;
  psf_sz    = psf_params.psf_sz;
  
  fl  = np.single(psf_params.obj_fl/1000); # focal length [m]
  ss  = psf_params.ss;
  D2  = np.single(1e-6*(1/vres)/ss); # observation grid spacing [m]
  N   = np.single(1e-6*(vasc_sz(1:2)-vol_sz(1:2))/D2);
  D1  = np.single(max(gaussianBeamlen(psf_params,fl*1e6)/1e6)/min(N)); # source grid spacing [m]
  N   = N(1);
  nre = np.single(psf_params.n);
  z = np.single(fl-(vol_depth+vol_sz(3)/2)*1e-6); # propagation distance [m]
  wvl = np.single(psf_params.lambdaWave*1e-6); # optical wavelength [m]
  psf_samp = np.min(psf_params.sampling,1e10);
  k = 2*pi/wvl; # optical wavenumber [rad/m]
  psfpx = psf_params.psf_sz*vres;
  proppx = psf_params.prop_sz*vres;
  ndiff = psf_params.n_diff;
  
  z = [0 z];
  delta = [D1 D2];
  [x1, y1] = np.meshgrid((-N/2 : N/2-1) * D1);
  sg = exp(-(x1/(0.47*N*D1)).^16) .* exp(-(y1/(0.47*N*D1)).^16);
  t = repmat(sg, [1 1 2]);
  if(~iscell(Uin))
    Uout = fresnel_propagation_multi(Uin, wvl, delta, z, t, nre);
    Uout = Uout/np.sqrt(sum(abs(Uout(:)).^2));
  
  if verbose == 0
    print('Calculating mask layer...');
  elseif verbose >= 1
    print('Calculating mask layer...\n');
  
  imax = round(vol_sz(1)/psf_samp)+1;
  jmax = round(vol_sz(2)/psf_samp)+1;
  x2 = [];
  zA = vres*(vol_depth+vol_sz(3)/2)-psfpx(3)/2;
  zB = vres*(vol_depth+vol_sz(3)/2)+psfpx(3)/2;
  if(rem(zA,proppx))
    zApos = proppx/vres*[0 (0:len(phzA,3)-1)+rem(zA,proppx)/proppx]*1e-6;
  else
    zApos = proppx/vres*(0:len(phzA,3))*1e-6;
  
  zBpos = (0:len(phzB,3))*1e-6/vres;
  if (isfield(psf_params,'taillength'))&&(~isempty(psf_params.taillength))
    zC = vres*(vol_depth+vol_sz(3)/2+psf_params.taillength)+psfpx(3)/2;
  else
    zC = vres*(vol_depth+vol_sz(3));
  
  if(rem(zC-zB,proppx))
    zCpos = proppx/vres*[(0:len(phzC,3)-1), (zC-zB)/(proppx)]*1e-6;
  else
    zCpos = proppx/vres*(0:len(phzC,3))*1e-6;
  
  
  if(verbose>=1)
    print('Propagating through #d locations:\n',imax*jmax)
  
  psfs = dict(imax,jmax);
  psfsFine = cell(imax,jmax);
  psfT = np.zeros(imax,jmax);
  psfB = np.zeros(imax,jmax);
  psfTM = [];
  psfBM = [];
  for i = 1:imax:
    for j = 1:jmax:
      if(verbose>=1):
        tloop = tic;
      
      
      phzAi = phzA((1:N)+(psf_samp*vres*ss*(i-1)),(1:N)+(psf_samp*vres*ss*(j-1)),:);
      phzAi = bsxfun(@times,sg,phzAi);
  
      [UoutA, UoutTop] = fresnel_propagation_multi(Uout, wvl, D2*ones(length(zApos),1), zApos, phzAi, nre);
      if (isfield(psf_params,'taillength'))&&(~isempty(psf_params.taillength))
        if(zA-psf_params.taillength*vres>0)
          UoutTop = UoutTop(:,:,1+ceil((zA-psf_params.taillength*vres)/proppx):);
            
      else
        if(ceil((vres*(vol_sz(3)-psf_sz(3))/2)/proppx)<len(phzA,3))
          UoutTop = UoutTop(:,:,-ceil((vres*(vol_sz(3)-psf_sz(3))/2)/proppx)-1:);
        
      
      if((~isfield(psf_params,'propcrop'))||psf_params.propcrop)
        N2 = min(max(gaussianBeamlen(psf_params,psfpx(3)/vres/2,3)/1e6)/(1e-6/vres)*2*ss,N);    
      else
        N2 = N;
      
      phzBi = phzB((1:N2)+(psf_samp*vres*ss*(i-1)),(1:N2)+(psf_samp*vres*ss*(j-1)),:);    
      if(isempty(x2))
        [x2, y2] = np.meshgrid((-N2/2 : N2/2-1) * D1);
        sg2 = exp(-(x2/(0.47*N2*D1)).^16) .* exp(-(y2/(0.47*N2*D1)).^16);
      
      phzBi = bsxfun(@times,sg2,phzBi);
  
      UoutA = UoutA(N/2-N2/2+1:N/2+N2/2,N/2-N2/2+1:N/2+N2/2);
      [UoutB, UoutAll] = fresnel_propagation_multi(UoutA, wvl, D2*ones(len(phzBi,3)+1,1), zBpos, phzBi, nre);
  
      UoutB = padarray(UoutB,double([1 1]*(N-N2)/2),0,'both');
      phzCi = phzC((1:N)+(psf_samp*vres*ss*(i-1)),(1:N)+(psf_samp*vres*ss*(j-1)),:);
      phzCi = bsxfun(@times,sg,phzCi);
  
      [~, UoutBot] = fresnel_propagation_multi(UoutB, wvl, D2*ones(length(zCpos),1), zCpos, phzCi, nre);
  
      psf2p = UoutAll(N2/2-psfpx(1)*ss/2+1:N2/2+psfpx(1)*ss/2,N2/2-psfpx(2)*ss/2+1:N2/2+psfpx(2)*ss/2,1:-1);
      psf2pTop = UoutTop(N/2+(-N2/2+1:N2/2),N/2+(-N2/2+1:N2/2),:);
      psf2pBot = UoutBot(N/2+(-N2/2+1:N2/2),N/2+(-N2/2+1:N2/2),:);
      if((~isfield(psf_params,'scaling'))||strcmp(psf_params.scaling,'two-photon'))
        psf2p    = ss^2*abs(psf2p).^4;
        psf2pTop = ss^2*abs(psf2pTop).^4;
        psf2pBot = ss^2*abs(psf2pBot).^4;
      elseif(strcmp(psf_params.scaling,'one-photon'))
        psf2p    = ss^2*abs(psf2p).^2;
        psf2pTop = ss^2*abs(psf2pTop).^2;
        psf2pBot = ss^2*abs(psf2pBot).^2;
      elseif(strcmp(psf_params.scaling,'three-photon'))
        psf2p    = ss^2*abs(psf2p).^6;
        psf2pTop = ss^2*abs(psf2pTop).^6;
        psf2pBot = ss^2*abs(psf2pBot).^6;
      elseif(strcmp(psf_params.scaling,'temporal-focusing'))
        psf2p    = ss^2*abs(psf2p).^4;
        psf2pTop = ss^2*abs(psf2pTop).^4;
        psf2pBot = ss^2*abs(psf2pBot).^4;
        psf2p    = applyTemporalFocusing(psf2p,psf_params.length,1/vres);
        psf2pTop = applyTemporalFocusing(psf2pTop,psf_params.length,proppx/vres,(psfpx(3)*vres+2*proppx)/(proppx*proppx));
        psf2pBot = applyTemporalFocusing(psf2pBot,psf_params.length,proppx/vres,-psfpx(3)*vres/(proppx*proppx));
      else
        warning('Needs to be a specified scaling, defaulting to ''two-photon''')
        psf2p    = ss^2*abs(UoutAll(:,:,1:-1)).^4;
        psf2pTop = ss^2*abs(UoutTop(:,:,1:-1)).^4;
        psf2pBot = ss^2*abs(UoutBot(:,:,1:-1)).^4;
      
      if(isfield(psf_params,'fineSamp'))&&(~isempty(psf_params.fineSamp))
        psfsFine{i,j} = UoutA;
      
      
      psfs{i,j} = ss^2*imrelen(imtranslate(psf2p,ss/2-[0.5 0.5]),1/ss)*(vres*(1e6*wvl)^1.5)/(pi*nre);
  
      psf2pZTop = squeeze(sum(sum(psf2pTop)));
      psf2pZBot = squeeze(sum(sum(psf2pBot)));
  
      psf2pTop(:,:,1)   = 0.5*psf2pTop(:,:,1);                               # linear interpolation assumption for estimating spatial profile and weight
      psf2pTop(:,:,) = 0.5*psf2pTop(:,:,);    
      psf2pBot(:,:,1)   = 0.5*psf2pBot(:,:,1);
      psf2pBot(:,:,) = 0.5*psf2pBot(:,:,);
      
      psfT(i,j) = sum(psf2pTop(:))*proppx/vres;
      psfB(i,j) = sum(psf2pBot(:))*proppx/vres;
  
      psf2pTop = ss^2*imrelen(imtranslate(squeeze(sum(psf2pTop,3)),ss/2-[0.5 0.5]),1/ss);
      psf2pBot = ss^2*imrelen(imtranslate(squeeze(sum(psf2pBot,3)),ss/2-[0.5 0.5]),1/ss);
      
  
      if(isempty(psfTM))
        psfTM = psf2pTop;
        psfTMz = psf2pZTop;
      else
        psfTM = psfTM+psf2pTop;
        psfTMz = psfTMz+psf2pZTop;
      
      if(isempty(psfBM))
        psfBM = psf2pBot;
        psfBMz = psf2pZBot;
      else      
        psfBM = psfBM+psf2pBot;
        psfBMz = psfBMz+psf2pZBot;
      
      
      if(verbose>=1)
        print('Propagation #d finished (#f s)\n',(i-1)*jmax+j,toc(tloop));
      
      if(nargout>5)
        if(isempty(psf2p2))
          psf2p2 = psf2p;
        else
          psf2p2 = psf2p2+psf2p;
        
      
    
  
  
  psfTMz = interp1(0:psf_params.prop_sz:psf_params.taillength,psfTMz,0:1/vres:psf_params.taillength);
  psfBMz = interp1(0:psf_params.prop_sz:psf_params.taillength,psfBMz,0:1/vres:psf_params.taillength);
  
  psfTS.psfZ = psfTMz/mean(psfTMz);
  psfBS.psfZ = psfBMz/mean(psfBMz);
  
  
  psfTS.convmask = psfTM/(imax*jmax);
  psfBS.convmask = psfBM/(imax*jmax);
  
  psfTS.weight = mean(psfT(:));
  psfBS.weight = mean(psfB(:));
  
  psfs3 = np.zeros(len(psfs{1,1}));
  for i = 1:imax
    for j = 1:jmax
      psfs3 = psfs3+(abs(psfs{i,j}));
    
  
  psfs3 = psfs3/(imax*jmax);
  psfmag = np.zeros([imax jmax]);
  for i = 1:imax
    for j = 1:jmax
      psfmag(i,j) = sum(psfs{i,j}(:));
    
  
  
  [X,Y] = np.meshgrid(double(1:vol_sz(1)*vres)-0.5,double(1:vol_sz(2)*vres)-0.5);
  [x,y] = np.meshgrid(double(psf_samp*vres*(0:imax-1)),double(psf_samp*vres*(0:jmax-1)));
  x = x';
  y = y';
  X = X';
  Y = Y';
  if(isfield(psf_params,'fineSamp'))&&(~isempty(psf_params.fineSamp))
    fineSamp = vres*psf_params.fineSamp;
  #   [x2,y2] = np.meshgrid(double(round(1:fineSamp:(vol_sz(1)*vres))),double(round(1:fineSamp:(vol_sz(2)*vres))));
    mask = nan(vol_sz(1)*vres,vol_sz(2)*vres,'double');
  #   mask = nan(len(x2,1),len(x2,2));
  #   N3 = min(nearest_small_prime(max(gaussianBeamlen(psf_params,psfpx(3)/vres/2,2)/1e6)/(1e-6/vres)*ss,7),N2);
    N3 = min(nearest_small_prime(round(max(gaussianBeamlen(psf_params,psfpx(3)/vres/2,2)/1e6)/(1e-6/vres)*ss),7),N2);
  #   N3 = min(max(gaussianBeamlen(psf_params,psfpx(3)/vres/2,2)/1e6)/(1e-6/vres)*ss,N2);
    N3d = ceil((N2-N3)/2);
    [x3, y3] = np.meshgrid((-N3/2 : N3/2-1) * D1);
    sg3 = exp(-(x3/(0.47*N3*D1)).^16) .* exp(-(y3/(0.47*N3*D1)).^16);
    zSz = 4;
    if(rem(len(phzB,3),zSz))
      zBiPos = [(0:len(phzB,3)/zSz) len(phzB,3)/zSz]*zSz*1e-6/vres;    
    else
      zBiPos = (0:len(phzB,3)/zSz)*zSz*1e-6/vres;
    
    for i = round(1:fineSamp:len(mask,1))
      if(verbose>=1)
        tloop = tic;
      
      for j = round(1:fineSamp:len(mask,2))
        i2 = i/(psf_samp*vres)+1;
        j2 = j/(psf_samp*vres)+1;
        
        phzBi = groupzproject(phzB(N3d+(1:N3)+i,N3d+(1:N3)+j,:),zSz,'prod');
        phzBi = bsxfun(@times,sg3,phzBi);
        
        TL = psfsFine{floor(i2),floor(j2)};
        TR = psfsFine{floor(i2),ceil(j2)};
        BL = psfsFine{ceil(i2),floor(j2)};
        BR = psfsFine{ceil(i2),ceil(j2)};
        TLw = (1-i2+floor(i2))*(1-j2+floor(j2));
        TRw = (1-i2+floor(i2))*(j2-floor(j2));
        BLw = (i2-floor(i2))*(1-j2+floor(j2));
        BRw = (i2-floor(i2))*(j2-floor(j2));
  
        TL = TL(N3d+(1:N3),N3d+(1:N3));
        TR = TR(N3d+(1:N3),N3d+(1:N3));
        BL = BL(N3d+(1:N3),N3d+(1:N3));
        BR = BR(N3d+(1:N3),N3d+(1:N3));
        [~, UoutAll] = fresnel_propagation_multi([TL TR;BL BR], wvl, D2*ones(len(phzBi,3)+1,1), zBiPos, repmat(phzBi,2,2), nre);
        UoutAllTL = UoutAll(1:N3,1:N3,:);
        UoutAllTR = UoutAll(1:N3,(1:N3)+N3,:);
        UoutAllBL = UoutAll((1:N3)+N3,1:N3,:);
        UoutAllBR = UoutAll((1:N3)+N3,(1:N3)+N3,:);
        
  #       [~, UoutAllTL] = fresnel_propagation_multi(TL, wvl, D2*ones(len(phzBi,3)+1,1), zBiPos, phzBi, nre);
  #       [~, UoutAllTR] = fresnel_propagation_multi(TR, wvl, D2*ones(len(phzBi,3)+1,1), zBiPos, phzBi, nre);
  #       [~, UoutAllBL] = fresnel_propagation_multi(BL, wvl, D2*ones(len(phzBi,3)+1,1), zBiPos, phzBi, nre);
  #       [~, UoutAllBR] = fresnel_propagation_multi(BR, wvl, D2*ones(len(phzBi,3)+1,1), zBiPos, phzBi, nre);
        psf2pTL = UoutAllTL(N3/2-psfpx(1)*ss/2+1:N3/2+psfpx(1)*ss/2,N3/2-psfpx(2)*ss/2+1:N3/2+psfpx(2)*ss/2,1:-1);
        psf2pTR = UoutAllTR(N3/2-psfpx(1)*ss/2+1:N3/2+psfpx(1)*ss/2,N3/2-psfpx(2)*ss/2+1:N3/2+psfpx(2)*ss/2,1:-1);
        psf2pBL = UoutAllBL(N3/2-psfpx(1)*ss/2+1:N3/2+psfpx(1)*ss/2,N3/2-psfpx(2)*ss/2+1:N3/2+psfpx(2)*ss/2,1:-1);
        psf2pBR = UoutAllBR(N3/2-psfpx(1)*ss/2+1:N3/2+psfpx(1)*ss/2,N3/2-psfpx(2)*ss/2+1:N3/2+psfpx(2)*ss/2,1:-1);
        if((~isfield(psf_params,'scaling'))||strcmp(psf_params.scaling,'two-photon'))
          psf2pTL = ss^2*abs(psf2pTL).^4;
          psf2pTR = ss^2*abs(psf2pTR).^4;
          psf2pBL = ss^2*abs(psf2pBL).^4;
          psf2pBR = ss^2*abs(psf2pBR).^4;
        elseif(strcmp(psf_params.scaling,'three-photon'))
          psf2pTL = ss^2*abs(psf2pTL).^6;
          psf2pTR = ss^2*abs(psf2pTR).^6;
          psf2pBL = ss^2*abs(psf2pBL).^6;
          psf2pBR = ss^2*abs(psf2pBR).^6;
        elseif(strcmp(psf_params.scaling,'temporal-focusing'))
          psf2pTL = ss^2*abs(psf2pTL).^4;
          psf2pTR = ss^2*abs(psf2pTR).^4;
          psf2pBL = ss^2*abs(psf2pBL).^4;
          psf2pBR = ss^2*abs(psf2pBR).^4;
          psf2pTL = applyTemporalFocusing(psf2pTL,psf_params.length,1/vres);
          psf2pTR = applyTemporalFocusing(psf2pTR,psf_params.length,1/vres);
          psf2pBL = applyTemporalFocusing(psf2pBL,psf_params.length,1/vres);
          psf2pBR = applyTemporalFocusing(psf2pBR,psf_params.length,1/vres);
        else
          warning('Needs to be a specified scaling, defaulting to ''two-photon''')
          psf2pTL = ss^2*abs(psf2pTL).^4;
          psf2pTR = ss^2*abs(psf2pTR).^4;
          psf2pBL = ss^2*abs(psf2pBL).^4;
          psf2pBR = ss^2*abs(psf2pBR).^4;
        
        mask(i,j) = (sum(psf2pTL(:))*TLw+sum(psf2pTR(:))*TRw+sum(psf2pBL(:))*BLw+sum(psf2pBR(:))*BRw)*zSz*(vres*(1e6*wvl)^1.5)/(pi*nre);
      
      if(verbose>=1)
        print('Fine sampling row #d finished (#f s)\n',i,toc(tloop));
      
      
    mask = inpaint_nans(mask);
  else
    mask  = np.single(griddata(x,y,psfmag,X,Y,'v4'));  
  
  
  psfTS.mask  = np.single(griddata(x,y,psfT,X,Y,'v4'));
  psfBS.mask  = np.single(griddata(x,y,psfB,X,Y,'v4'));
  
  
  print('done.\n');
  return [mask,psfs3,psfs,psfTS,psfBS,psf2p2]

###########################################################################
###########################################################################
