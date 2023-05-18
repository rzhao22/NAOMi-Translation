import numpy as np

def constrainEstToSomas(est,sL):

# constrainEstToSomas(est,sL)
# return est2
# Function to constrain the est class to a subset of profiles
# that represent, for example, somatic components.
# 
# 2020 - Adam Charles

    tmpSL     = np.full((est.estactIdxs.shape), False)
    tmpSL[sL] = True
########################################################################################################################
    ## First: indexing hell
    class Est2:
        pass
    est2 = Est2()
    est2.estactIdxs      = tmpSL&est.estactIdxs                             # Full logical list of global components to keep
    est2.estactidealIdxs = tmpSL&est.estactidealIdxs                        # |- Same for ideal 

    est2.Idxs            = np.empty(est2.estactIdxs.shape).fill(np.nan)     # Create nan list
    est2.idealIdxs       = np.empty(est2.estactidealIdxs.shape).fill(np.nan)# |- Same for ideal
    est2.Idxs[est2.estactIdxs]           = range(sum(est2.estactIdxs))      # Populate with index of component in new list (for pairing purposes
    est2.idealIdxs[est2.estactidealIdxs] = range(sum(est2.estactidealIdxs)) # |- Same for ideal

    ## Need indexing to subselect new components from previously wittled list

    keepIdxs       = est2.Idxs[est.estactIdxs]
    keepIdealIdxs  = est2.idealIdxs[est.estactidealIdxs]

    est2.estact        = np.concatenate((est.estact[np.invert(np.isnan(keepIdxs)),:], est.estact[-1,:]), axis = 0)
    est2.estactideal   = np.concatenate((est.estactideal[np.invert(np.isnan(keepIdealIdxs)),:], est.estactideal[-1,:]), axis = 0)
    est2.compsIdealAll = est.compsIdealAll
    est2.compsIdeal    = est.compsIdeal[:,:,np.invert(np.isnan(keepIdxs))]

    est2.corrIdxs      = est.corrIdxs
    est2.corrvals      = est.corrvals
    est2.corrvalsIdeal = est.corrvalsIdeal
    est2.corrvalsIdeal[np.invert(np.isnan(est2.idealIdxs))] = float('nan')
    est2.corrSortIdeal = np.sort(est2.corrvalsIdeal,'descend')
    IXC                = np.argsort(est2.corrvalsIdeal,'descend')
    est2.corrIdxsIdeal       = est2.idealIdxs[IXC]
    return est2

# x - estact
# x - estactIdxs
# x - estactideal
# x - estactidealIdxs
# x - Idxs
# x - idealIdxs
# x - corrvals
# x - corrvalsIdeal
# x - corrSort
# x - corrIdxs
# x - corrSortIdeal
# x - corrIdxsIdeal
# x - compsIdealAll
# x - compsIdeal
# o - pairs 
# o - simTraces
# o - pairedTraces
# o - allpairs
# o - strongbound
# o - weakbound
# o - upbound
