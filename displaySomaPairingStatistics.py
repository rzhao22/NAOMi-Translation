def displayPairingStatistics(cnmf,pcaica,suite2p,est,varargin)


    if nargin > 4
        somaList = varargin{1}

    sL = false(size(est.corrvals))
    sL(somaList) = true

    size(sL)
    size(est.corrvals)

    estCV1 = (est.corrvals>0.1)&sL
    estCV3 = (est.corrvals>0.3)&sz
    estCV5 = (est.corrvals>0.5)&sL

    estCVP1 = unique(est.pairs(estCV1,1))
    estCVP3 = unique(est.pairs(estCV3,1))
    estCVP5 = unique(est.pairs(estCV5,1))

    cnmfCV1 = (cnmf.corrvals>0.1)&sL
    cnmfCV3 = (cnmf.corrvals>0.3)&sL
    cnmfCV5 = (cnmf.corrvals>0.5)&sL

    cnmfCVP1 = unique(cnmf.pairs(cnmfCV1,1))
    cnmfCVP3 = unique(cnmf.pairs(cnmfCV3,1))
    cnmfCVP5 = unique(cnmf.pairs(cnmfCV5,1))

    suite2pCV1 = (suite2p.corrvals>0.1)&sL
    suite2pCV3 = (suite2p.corrvals>0.3)&sL
    suite2pCV5 = (suite2p.corrvals>0.5)&sL

    suite2pCVP1 = unique(suite2p.pairs(suite2pCV1,1))
    suite2pCVP3 = unique(suite2p.pairs(suite2pCV3,1))
    suite2pCVP5 = unique(suite2p.pairs(suite2pCV5,1))

    pcaicaCV1 = (pcaica.corrvals>0.1)&sL
    pcaicaCV3 = (pcaica.corrvals>0.3)&sL
    pcaicaCV5 = (pcaica.corrvals>0.5)&sL

    pcaicaCVP1 = unique(pcaica.pairs(pcaicaCV1,1))
    pcaicaCVP3 = unique(pcaica.pairs(pcaicaCV3,1))
    pcaicaCVP5 = unique(pcaica.pairs(pcaicaCV5,1))


    print('Ideal:')
    print('Number componants: %d\n', size(est.estact,1))
    print('Number paired componants (>0.1,>0.3,>0.5): (%d,%d,%d)\n', sum(estCV1),sum(estCV3),sum(estCV5))
    print('Number unique componants found (>0.1,>0.3,>0.5): (%d,%d,%d)\n', 
        numel(estCVP1,1), numel(estCVP3,1), numel(estCVP5,1))
    print('Number doubled componants (>0.1,>0.3,>0.5): (%d,%d,%d)\n', sum((hist(est.pairs(estCV1,1),estCVP1))>1),
        sum((hist(est.pairs(estCV3,1),estCVP3))>1), sum((hist(est.pairs(estCV5,1),estCVP5))>1))

    print('CNMF:')
    print('Number componants: %d\n', size(cnmf.compFluoresence,1))
    print('Number paired componants (>0.1,>0.3,>0.5): (%d,%d,%d)\n', sum(cnmfCV1),sum(cnmfCV3),sum(cnmfCV5))
    print('Number unique componants found (>0.1,>0.3,>0.5): (%d,%d,%d)\n', 
        numel(cnmfCVP1), numel(cnmfCVP3), numel(cnmfCVP5))
    print('Number doubled componants (>0.1,>0.3,>0.5): (%d,%d,%d)\n', 
        sum((hist(cnmf.pairs(cnmfCV1,1),cnmfCVP1))>1),
        sum((hist(cnmf.pairs(cnmfCV3,1),cnmfCVP3))>1),
        sum((hist(cnmf.pairs(cnmfCV5,1),cnmfCVP5))>1))

    print('Suite2p:')
    print('Number componants: %d\n', size(suite2p.compTimecourse,1))
    print('Number paired componants (>0.1,>0.3,>0.5): (%d,%d,%d)\n', sum(suite2pCV1),sum(suite2pCV3),sum(suite2pCV5))
    print('Number unique componants found (>0.1,>0.3,>0.5): (%d,%d,%d)\n', 
        numel(suite2pCVP1), numel(suite2pCVP3), numel(suite2pCVP5))
    print('Number doubled componants (>0.1,>0.3,>0.5): (%d,%d,%d)\n', 
        sum((hist(suite2p.pairs(suite2pCV1,1),suite2pCVP1))>1),
        sum((hist(suite2p.pairs(suite2pCV3,1),suite2pCVP3))>1),
        sum((hist(suite2p.pairs(suite2pCV5,1),suite2pCVP5))>1))

    print('PCA/ICA:')
    print('Number componants: %d\n', size(pcaica.compTimecourse,1))
    print('Number paired componants (>0.1,>0.3,>0.5): (%d,%d,%d)\n', sum(pcaicaCV1),sum(pcaicaCV3),sum(pcaicaCV5))
    print('Number unique componants found (>0.1,>0.3,>0.5): (%d,%d,%d)\n', 
        numel(pcaicaCVP1), numel(pcaicaCVP3), numel(pcaicaCVP5))
    print('Number doubled componants (>0.1,>0.3,>0.5): (%d,%d,%d)\n', 
        sum((hist(pcaica.pairs(pcaicaCV1,1),pcaicaCVP1))>1),
        sum((hist(pcaica.pairs(pcaicaCV3,1),pcaicaCVP3))>1),
        sum((hist(pcaica.pairs(pcaicaCV5,1),pcaicaCVP5))>1))

    TMPcs1 = intersect(cnmfCVP1,   suite2pCVP1)
    TMPcp1 = intersect(cnmfCVP1,   pcaicaCVP1 )
    TMPsp1 = intersect(pcaicaCVP1, suite2pCVP1)
    TMPxx1 = intersect(TMPcs1,intersect(TMPcp1,TMPsp1))

    TMPcs3 = intersect(cnmfCVP3,   suite2pCVP3)
    TMPcp3 = intersect(cnmfCVP3,   pcaicaCVP3 )
    TMPsp3 = intersect(pcaicaCVP3, suite2pCVP3)
    TMPxx3 = intersect(TMPcs3,intersect(TMPcp3,TMPsp3))

    TMPcs5 = intersect(cnmfCVP5,   suite2pCVP5)
    TMPcp5 = intersect(cnmfCVP5,   pcaicaCVP5 )
    TMPsp5 = intersect(pcaicaCVP5, suite2pCVP5)
    TMPxx5 = intersect(TMPcs5,intersect(TMPcp5,TMPsp5))

    print('Number of unique cells found by:\n')
    print('\t - Only CNMF (>0.1,>0.3,>0.5): %d, %d, %d\n', 
        numel(setdiff(cnmfCVP1,    union(TMPcs1,TMPcp1))), 
        numel(setdiff(cnmfCVP3,    union(TMPcs3,TMPcp3))), 
        numel(setdiff(cnmfCVP5,    union(TMPcs5,TMPcp5))))
    print('\t - Only Suite2p (>0.1,>0.3,>0.5): %d, %d, %d\n', 
        numel(setdiff(suite2pCVP1, union(TMPcs1,TMPsp1))), 
        numel(setdiff(suite2pCVP3, union(TMPcs3,TMPsp3))), 
        numel(setdiff(suite2pCVP5, union(TMPcs5,TMPsp5))))
    print('\t - Only PCAICA (>0.1,>0.3,>0.5): %d, %d, %d\n', 
        numel(setdiff(pcaicaCVP1,  union(TMPcs1,TMPsp1))), 
        numel(setdiff(pcaicaCVP3,  union(TMPcs3,TMPsp3))), 
        numel(setdiff(pcaicaCVP5,  union(TMPcs5,TMPsp5))))
    print('\t - CNMF and Suite2p (>0.1,>0.3,>0.5): %d, %d, %d\n',    numel(TMPcs1), numel(TMPcs3), numel(TMPcs5))
    print('\t - CNMF and PCA/ICA (>0.1,>0.3,>0.5): %d, %d, %d\n',    numel(TMPcp1), numel(TMPcp3), numel(TMPcp5))
    print('\t - PCA/ICA and Suite2p (>0.1,>0.3,>0.5): %d, %d, %d\n', numel(TMPsp1), numel(TMPsp3), numel(TMPsp5))
    print('\t - All algorithms (>0.1,>0.3,>0.5): %d, %d, %d\n',      numel(TMPxx1), numel(TMPxx3), numel(TMPxx5))


    TMPinc1 = setdiff(estCVP1, cnmfCVP1   )
    TMPins1 = setdiff(estCVP1, suite2pCVP1)
    TMPinp1 = setdiff(estCVP1, pcaicaCVP1 )
    TMPinc3 = setdiff(estCVP3, cnmfCVP3   )
    TMPins3 = setdiff(estCVP3, suite2pCVP3)
    TMPinp3 = setdiff(estCVP3, pcaicaCVP3 )
    TMPinc5 = setdiff(estCVP5, cnmfCVP5   )
    TMPins5 = setdiff(estCVP5, suite2pCVP5)
    TMPinp5 = setdiff(estCVP5, pcaicaCVP5 )

    print('Number of unique cells found by the Ideal components:\n')
    print('\t - and NOT CNMF: (>0.1,>0.3,>0.5): %d, %d, %d\n',    numel(TMPinc1), numel(TMPinc3), numel(TMPinc5))
    print('\t - and NOT Suite2p: (>0.1,>0.3,>0.5): %d, %d, %d\n', numel(TMPins1), numel(TMPins3), numel(TMPins5))
    print('\t - and NOT PCA/ICA: (>0.1,>0.3,>0.5): %d, %d, %d\n', numel(TMPinp1), numel(TMPinp3), numel(TMPinp5))

    TMPcni1 = setdiff(cnmfCVP1,    estCVP1)
    TMPsni1 = setdiff(suite2pCVP1, estCVP1)
    TMPpni1 = setdiff(pcaicaCVP1,  estCVP1)
    TMPcni3 = setdiff(cnmfCVP3,    estCVP3)
    TMPsni3 = setdiff(suite2pCVP3, estCVP3)
    TMPpni3 = setdiff(pcaicaCVP3,  estCVP3)
    TMPcni5 = setdiff(cnmfCVP5,    estCVP5)
    TMPsni5 = setdiff(suite2pCVP5, estCVP5)
    TMPpni5 = setdiff(pcaicaCVP5,  estCVP5)

    print('Number of unique cells not found by the Ideal components:\n')
    print('\t - but found by CNMF: (>0.1,>0.3,>0.5): %d, %d, %d\n',    numel(TMPcni1), numel(TMPcni3), numel(TMPcni5))
    print('\t - but found by Suite2p: (>0.1,>0.3,>0.5): %d, %d, %d\n', numel(TMPsni1), numel(TMPsni3), numel(TMPsni5))
    print('\t - but found by PCA/ICA: (>0.1,>0.3,>0.5): %d, %d, %d\n', numel(TMPpni1), numel(TMPpni3), numel(TMPpni5))


