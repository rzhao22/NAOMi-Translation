import numpy as np
import matplotlib.pyplot as plt


# function displayPairingStatistics(cnmf,pcaica,suite2p,est,varargin)
def displayPairingStatistics(cnmf, pcaica, suite2p, est, *args):


    if len(args) > 0:
        somaList = args[0]

    sL           = np.full((est.corrvals.shape), False)
    sL[somaList] = True

    #clc

    print(sL.shape)
    print(est.corrvals.shape)

    estCV1 = (est.corrvals>0.1)&sL
    estCV3 = (est.corrvals>0.3)&sL
    estCV5 = (est.corrvals>0.5)&sL

    estCVP1 = np.unique(est.pairs[estCV1,0])
    estCVP3 = np.unique(est.pairs[estCV3,0])
    estCVP5 = np.unique(est.pairs[estCV5,0])

    cnmfCV1 = (cnmf.corrvals>0.1)&sL
    cnmfCV3 = (cnmf.corrvals>0.3)&sL
    cnmfCV5 = (cnmf.corrvals>0.5)&sL

    cnmfCVP1 = np.unique(cnmf.pairs[cnmfCV1,0])
    cnmfCVP3 = np.unique(cnmf.pairs[cnmfCV3,0])
    cnmfCVP5 = np.unique(cnmf.pairs[cnmfCV5,0])

    suite2pCV1 = (suite2p.corrvals>0.1)&sL
    suite2pCV3 = (suite2p.corrvals>0.3)&sL
    suite2pCV5 = (suite2p.corrvals>0.5)&sL

    suite2pCVP1 = np.unique(suite2p.pairs[suite2pCV1,0])
    suite2pCVP3 = np.unique(suite2p.pairs[suite2pCV3,0]) 
    suite2pCVP5 = np.unique(suite2p.pairs[suite2pCV5,0]) 

    pcaicaCV1 = (pcaica.corrvals>0.1)&sL
    pcaicaCV3 = (pcaica.corrvals>0.3)&sL
    pcaicaCV5 = (pcaica.corrvals>0.5)&sL

    pcaicaCVP1 = np.unique(pcaica.pairs[pcaicaCV1,0])
    pcaicaCVP3 = np.unique(pcaica.pairs[pcaicaCV3,0])
    pcaicaCVP5 = np.unique(pcaica.pairs[pcaicaCV5,0])

##########################################################################%
## Single algorithm statistics

    print('Ideal:')
    print('Number componants: ', np.size(est.estact,0))
    print('Number paired componants (>0.1,>0.3,>0.5): (', sum(estCV1),',',sum(estCV3),',',sum(estCV5),')')
    ##### numel(sth, 1) = 1?????????????????????????????????????????????????????????????????????????????????
    print('Number np.unique componants found (>0.1,>0.3,>0.5): (',
                         np.size(estCVP1,0), ',',np.size(estCVP3,0), ',',np.size(estCVP5,0),')')
    print('Number doubled componants (>0.1,>0.3,>0.5): (',
        sum((plt.hist(est.pairs(estCV1,0),estCVP1)[0])>1),',',
        sum((plt.hist(est.pairs(estCV3,0),estCVP3)[0])>1),',',
        sum((plt.hist(est.pairs(estCV5,0),estCVP5)[0])>1),')')

    print('CNMF:')
    print('Number componants: ', np.size(cnmf.compFluoresence,0))
    print('Number paired componants (>0.1,>0.3,>0.5): (', sum(cnmfCV1),',',sum(cnmfCV3),',',sum(cnmfCV5),')')
    print('Number np.unique componants found (>0.1,>0.3,>0.5): (',
                             np.size(cnmfCVP1), ',',np.size(cnmfCVP3), ',',np.size(cnmfCVP5),')')
    print('Number doubled componants (>0.1,>0.3,>0.5): (',
        sum((plt.hist(cnmf.pairs(cnmfCV1,0),cnmfCVP1)[0])>1),',',
        sum((plt.hist(cnmf.pairs(cnmfCV3,0),cnmfCVP3)[0])>1),',',
        sum((plt.hist(cnmf.pairs(cnmfCV5,0),cnmfCVP5)[0])>1),')')

    print('Suite2p:')
    print('Number componants: ', np.size(suite2p.compTimecourse,0))
    print('Number paired componants (>0.1,>0.3,>0.5): (', sum(suite2pCV1),',',sum(suite2pCV3),',',sum(suite2pCV5),')')
    print('Number np.unique componants found (>0.1,>0.3,>0.5): (',
                             np.size(suite2pCVP1), ',',np.size(suite2pCVP3), ',',np.size(suite2pCVP5),')')
    print('Number doubled componants (>0.1,>0.3,>0.5): (',
        sum((plt.hist(suite2p.pairs(suite2pCV1,0),suite2pCVP1)[0])>1),',',
        sum((plt.hist(suite2p.pairs(suite2pCV3,0),suite2pCVP3)[0])>1),',',
        sum((plt.hist(suite2p.pairs(suite2pCV5,0),suite2pCVP5)[0])>1),')')

    print('PCA/ICA:')
    print('Number componants: ', np.size(pcaica.compTimecourse,0))
    print('Number paired componants (>0.1,>0.3,>0.5): (', sum(pcaicaCV1),',',sum(pcaicaCV3),',',sum(pcaicaCV5),')')
    print('Number np.unique componants found (>0.1,>0.3,>0.5): (',
                             np.size(pcaicaCVP1), ',',np.size(pcaicaCVP3), ',',np.size(pcaicaCVP5),')')
    print('Number doubled componants (>0.1,>0.3,>0.5): (',
        sum((plt.hist(pcaica.pairs(pcaicaCV1,0),pcaicaCVP1))>1),',',
        sum((plt.hist(pcaica.pairs(pcaicaCV3,0),pcaicaCVP3))>1),',',
        sum((plt.hist(pcaica.pairs(pcaicaCV5,0),pcaicaCVP5))>1),')')


##########################################################################%
## Comparison statistics

    TMPcs1 = np.intersect1d(cnmfCVP1,   suite2pCVP1)
    TMPcp1 = np.intersect1d(cnmfCVP1,   pcaicaCVP1 )
    TMPsp1 = np.intersect1d(pcaicaCVP1, suite2pCVP1)
    TMPxx1 = np.intersect1d(TMPcs1,np.intersect1d(TMPcp1,TMPsp1))

    TMPcs3 = np.intersect1d(cnmfCVP3,   suite2pCVP3)
    TMPcp3 = np.intersect1d(cnmfCVP3,   pcaicaCVP3 )
    TMPsp3 = np.intersect1d(pcaicaCVP3, suite2pCVP3)
    TMPxx3 = np.intersect1d(TMPcs3,np.intersect1d(TMPcp3,TMPsp3))

    TMPcs5 = np.intersect1d(cnmfCVP5,   suite2pCVP5)
    TMPcp5 = np.intersect1d(cnmfCVP5,   pcaicaCVP5 )
    TMPsp5 = np.intersect1d(pcaicaCVP5, suite2pCVP5)
    TMPxx5 = np.intersect1d(TMPcs5,np.intersect1d(TMPcp5,TMPsp5))

    #TMPsp1 = np.intersect1d(pcaicaCVP5, suite2pCVP5)
    #TMPxxI = np.intersect1d(TMPcs5,np.intersect1d(TMPcp5,TMPsp5))

    print('Number of np.unique cells found by:')
    print('\t - Only CNMF (>0.1,>0.3,>0.5): ',
        np.size(np.setdiff1d(cnmfCVP1,    union(TMPcs1,TMPcp1))),',',
        np.size(np.setdiff1d(cnmfCVP3,    union(TMPcs3,TMPcp3))),',',
        np.size(np.setdiff1d(cnmfCVP5,    union(TMPcs5,TMPcp5))))
    print('\t - Only Suite2p (>0.1,>0.3,>0.5): ',
        np.size(np.setdiff1d(suite2pCVP1, union(TMPcs1,TMPsp1))),',',
        np.size(np.setdiff1d(suite2pCVP3, union(TMPcs3,TMPsp3))),',',
        np.size(np.setdiff1d(suite2pCVP5, union(TMPcs5,TMPsp5))))
    print('\t - Only PCAICA (>0.1,>0.3,>0.5): ',
        np.size(np.setdiff1d(pcaicaCVP1,  union(TMPcs1,TMPsp1))),',',
        np.size(np.setdiff1d(pcaicaCVP3,  union(TMPcs3,TMPsp3))),',',
        np.size(np.setdiff1d(pcaicaCVP5,  union(TMPcs5,TMPsp5))))
    print('\t - CNMF and Suite2p (>0.1,>0.3,>0.5): ',    np.size(TMPcs1),',', np.size(TMPcs3),',', np.size(TMPcs5))
    print('\t - CNMF and PCA/ICA (>0.1,>0.3,>0.5): ',    np.size(TMPcp1),',', np.size(TMPcp3),',', np.size(TMPcp5))
    print('\t - PCA/ICA and Suite2p (>0.1,>0.3,>0.5): ', np.size(TMPsp1),',', np.size(TMPsp3),',', np.size(TMPsp5))
    print('\t - All algorithms (>0.1,>0.3,>0.5): ',      np.size(TMPxx1),',', np.size(TMPxx3),',', np.size(TMPxx5))

    ## Comparison to ideal

    # TMPic1 = np.intersect1d(np.unique(est.pairs(est.corrvals>0.1,1)), np.unique(suite2p.pairs(suite2p.corrvals>0.1,1)));
    # TMPis1 = np.intersect1d(np.unique(est.pairs(est.corrvals>0.1,1)), np.unique(suite2p.pairs(suite2p.corrvals>0.1,1)));
    # TMPip1 = np.intersect1d(np.unique(est.pairs(est.corrvals>0.1,1)), np.unique(pcaica.pairs(pcaica.corrvals  >0.1,1)));
    # TMPic3 = np.intersect1d(np.unique(est.pairs(est.corrvals>0.3,1)), np.unique(suite2p.pairs(suite2p.corrvals>0.3,1)));
    # TMPis3 = np.intersect1d(np.unique(est.pairs(est.corrvals>0.3,1)), np.unique(suite2p.pairs(suite2p.corrvals>0.3,1)));
    # TMPip3 = np.intersect1d(np.unique(est.pairs(est.corrvals>0.3,1)), np.unique(pcaica.pairs(pcaica.corrvals  >0.3,1)));
    # TMPic5 = np.intersect1d(np.unique(est.pairs(est.corrvals>0.5,1)), np.unique(suite2p.pairs(suite2p.corrvals>0.5,1)));
    # TMPis5 = np.intersect1d(np.unique(est.pairs(est.corrvals>0.5,1)), np.unique(suite2p.pairs(suite2p.corrvals>0.5,1)));
    # TMPip5 = np.intersect1d(np.unique(est.pairs(est.corrvals>0.5,1)), np.unique(pcaica.pairs(pcaica.corrvals  >0.5,1)));

    TMPinc1 = np.setdiff1d(estCVP1, cnmfCVP1   )
    TMPins1 = np.setdiff1d(estCVP1, suite2pCVP1)
    TMPinp1 = np.setdiff1d(estCVP1, pcaicaCVP1 )
    TMPinc3 = np.setdiff1d(estCVP3, cnmfCVP3   )
    TMPins3 = np.setdiff1d(estCVP3, suite2pCVP3)
    TMPinp3 = np.setdiff1d(estCVP3, pcaicaCVP3 )
    TMPinc5 = np.setdiff1d(estCVP5, cnmfCVP5   )
    TMPins5 = np.setdiff1d(estCVP5, suite2pCVP5)
    TMPinp5 = np.setdiff1d(estCVP5, pcaicaCVP5 )

    print('Number of np.unique cells found by the Ideal components:')
    print('\t - and NOT CNMF: (>0.1,>0.3,>0.5): ',    np.size(TMPinc1),',', np.size(TMPinc3),',', np.size(TMPinc5))
    print('\t - and NOT Suite2p: (>0.1,>0.3,>0.5): ', np.size(TMPins1),',', np.size(TMPins3),',', np.size(TMPins5))
    print('\t - and NOT PCA/ICA: (>0.1,>0.3,>0.5): ', np.size(TMPinp1),',', np.size(TMPinp3),',', np.size(TMPinp5))

    TMPcni1 = np.setdiff1d(cnmfCVP1,    estCVP1)
    TMPsni1 = np.setdiff1d(suite2pCVP1, estCVP1)
    TMPpni1 = np.setdiff1d(pcaicaCVP1,  estCVP1)
    TMPcni3 = np.setdiff1d(cnmfCVP3,    estCVP3)
    TMPsni3 = np.setdiff1d(suite2pCVP3, estCVP3)
    TMPpni3 = np.setdiff1d(pcaicaCVP3,  estCVP3)
    TMPcni5 = np.setdiff1d(cnmfCVP5,    estCVP5)
    TMPsni5 = np.setdiff1d(suite2pCVP5, estCVP5)
    TMPpni5 = np.setdiff1d(pcaicaCVP5,  estCVP5)

    print('Number of np.unique cells not found by the Ideal components:')
    print('\t - but found by CNMF: (>0.1,>0.3,>0.5): ',    np.size(TMPcni1),',', np.size(TMPcni3),',', np.size(TMPcni5))
    print('\t - but found by Suite2p: (>0.1,>0.3,>0.5): ', np.size(TMPsni1),',', np.size(TMPsni3),',', np.size(TMPsni5))
    print('\t - but found by PCA/ICA: (>0.1,>0.3,>0.5): ', np.size(TMPpni1),',', np.size(TMPpni3),',', np.size(TMPpni5))


###########################################################################
###########################################################################

def union(lst1, lst2):
    final_list = list(set(lst1) | set(lst2))
    return final_list