import numpy as np
import vec
# genAlgorithmPairROC(cnmf,pcaica,suite2p,est,*args)
# return roc
def genAlgorithmPairROC(cnmf, pcaica, suite2p, est, *args):
    ###########################################################################
    ## Input parsing

    nargin = len(args)
    if nargin > 0:
        N = args[0]
    else:
        N = 0.5

    if nargin > 1:
        cutoff = args[1]
    else:
        cutoff = 0.5

    if not cutoff:
        cutoff = 0.5

    if not N:
        N = 100

    ###########################################################################
    ## Initializations

    ord_thresh = np.linspace(0,1,N)

    cnmfTT    = vec(sum(sum(cnmf.compSpatial))) * cnmf.compTimecourse
    pcaicaTT  = vec(sum(sum(pcaica.compSpatial))) * pcaica.compTimecourse
    estTT     = vec(sum(sum(est.compsIdeal))) * est.estact[:-2,:]
    suite2pTT = vec(sum(sum(suite2p.compSpatial > 0))) * (-np.median(suite2p.compTimecourse,axis = 1) + suite2p.compTimecourse)
    # suite2pTT = bsxfun(@times, vec(sum(sum(suite2p.compSpatial>0))), suite2p.compTimecourse);
    # suite2pTT = bsxfun(@plus,  -median(suite2pTT,2), suite2pTT);

    # define an empty class for storing roc
    class Roc:
        pass
    roc = Roc()

    roc.max.cnmf      = max(max(cnmfTT))
    roc.max.suite2p   = max(max(suite2pTT))
    roc.max.pcaica    = max(max(pcaicaTT))
    roc.max.est       = max(max(estTT))

    roc.min.cnmf      = min(np.amax(cnmfTT,axis = 1))
    roc.min.suite2p   = min(np.amax(suite2pTT,axis = 1))
    roc.min.pcaica    = min(np.amax(pcaicaTT,axis = 1))
    roc.min.est       = min(np.amax(estTT,axis = 1))


    roc.Ntrue.cnmf    = sum(cnmf.corrvals>cutoff)
    roc.Ntrue.suite2p = sum(suite2p.corrvals>cutoff)
    roc.Ntrue.pcaica  = sum(pcaica.corrvals>cutoff)
    roc.Ntrue.est     = sum(est.corrvals>cutoff)

    roc.Nflse.cnmf    = np.size(cnmfTT,0) - sum(cnmf.corrvals>cutoff)
    roc.Nflse.suite2p = np.size(suite2pTT,0) - sum(suite2p.corrvals>cutoff)
    roc.Nflse.pcaica  = np.size(pcaicaTT,0) - sum(pcaica.corrvals>cutoff)
    roc.Nflse.est     = np.size(estTT,0) - sum(est.corrvals>cutoff)

    roc.FA.cnmf    = np.zeros([np.size(ord_thresh),1])
    roc.FA.suite2p = np.zeros([np.size(ord_thresh),1])
    roc.FA.pcaica  = np.zeros([np.size(ord_thresh),1])
    roc.FA.est     = np.zeros([np.size(ord_thresh),1])

    roc.TP.cnmf    = np.zeros([np.size(ord_thresh),1])
    roc.TP.suite2p = np.zeros([np.size(ord_thresh),1])
    roc.TP.pcaica  = np.zeros([np.size(ord_thresh),1])
    roc.TP.est     = np.zeros([np.size(ord_thresh),1])

    ###########################################################################
    ## calculate ROC values

    for ll in range(np.size(ord_thresh)):
        # Pick time traces above a certain threshold
        roc.pick.cnmf    = np.where(np.amax(cnmfTT,axis = 1)    >= roc.min.cnmf    + ord_thresh(ll)*(roc.max.cnmf-roc.min.cnmf))
        roc.pick.suite2p = np.where(np.amax(suite2pTT,axis = 1) >= roc.min.suite2p + ord_thresh(ll)*(roc.max.suite2p-roc.min.suite2p))
        roc.pick.pcaica  = np.where(np.amax(pcaicaTT,axis = 1)  >= roc.min.pcaica  + ord_thresh(ll)*(roc.max.pcaica-roc.min.pcaica))
        roc.pick.est     = np.where(np.amax(estTT,axis = 1)     >= roc.min.est     + ord_thresh(ll)*(roc.max.est-roc.min.est))
        
        cnmf_HIT    = np.intersect1d(roc.pick.cnmf,    cnmf.pairs(cnmf.corrvals>cutoff,1))
        suite2p_HIT = np.intersect1d(roc.pick.suite2p, suite2p.pairs(suite2p.corrvals>cutoff,1))
        pcaica_HIT  = np.intersect1d(roc.pick.pcaica,  pcaica.pairs(pcaica.corrvals>cutoff,1))
        est_HIT     = np.intersect1d(roc.pick.est,     est.pairs(est.corrvals>cutoff,1))
        
        roc.TP.cnmf[ll]    = np.size(cnmf_HIT)/roc.Ntrue.cnmf
        roc.TP.suite2p[ll] = np.size(suite2p_HIT)/roc.Ntrue.suite2p
        roc.TP.pcaica[ll]  = np.size(pcaica_HIT)/roc.Ntrue.pcaica
        roc.TP.est[ll]     = np.size(est_HIT)/roc.Ntrue.est
        
        roc.FA.cnmf[ll]    = np.size(np.setdiff1d(roc.pick.cnmf,cnmf_HIT))/roc.Nflse.cnmf
        roc.FA.suite2p[ll] = np.size(np.setdiff1d(roc.pick.suite2p,suite2p_HIT))/roc.Nflse.suite2p
        roc.FA.pcaica[ll]  = np.size(np.setdiff1d(roc.pick.pcaica,pcaica_HIT))/roc.Nflse.pcaica
        roc.FA.est[ll]     = np.size(np.setdiff1d(roc.pick.est,est_HIT))/roc.Nflse.est

    return roc

###########################################################################
###########################################################################