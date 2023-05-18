
def correlateAllSegmentations(cnmf, pcaica, suite2p, est, neur_act2, *args):

    if len(args)!=0:
        sL = args[0]
    else:
        sL = np.arange(0,len(neur_act2))

    if len(args)>1:
        calc_only = args[1]
    else:
        calc_only = [True, True, True, True]



    if calc_only[0]:
        # calculate pairs based on overlap and best temporal correlation
        cnmf['corrvals'], cnmf['pairs'], cnmf['corrM'] = corrTimeAndSpace((neur_act2[sL]).T,
                                                                          (cnmf['compFlourescence']).T,
                                                                          '''est['compsIdealAll']''',
                                                                          cnmf['compSpatial'])
        # extract the paired traces (simulated) 
        cnmf['simTraces']= neur_act2[cnmf['pairs'][:,0],:]
        
        # extract the paired traces (estimated)
        cnmf['pairedTraces'] = cnmf['compFlouresence'][cnmf['pairs'][:,1],:]

    if calc_only[1]:
        # calculate pairs based on overlap and best temporal correlation
        pcaica['corrvals'], pcaica['pairs'], pcaica['corrM'] = corrTimeAndSpace((neur_act2[sL]).T,
                                                                                pcaica['compTimecourse'].T,
                                                                                '''est['compsIdealAll']''',
                                                                                pcaica['compSpatialSc'])                                                  
        # extract the paired traces (simulated)                                                                 
        pcaica['simTraces'] = neur_act2[pcaica['pairs'][:,0],:]

        # extract the paired traces (estimated)
        pcaica['pairedTraces'] = pcaica['compTimecourse'][pcaica['pairs'][:,1],:]

    if calc_only[2]:
        # calculate pairs based on overlap and best temporal correlation
        suite2p['corrvals'], suite2p['pairs'], suite2p['corrM'] = corrTimeAndSpace((neur_act2[sL]).T,
                                                                                suite2p['compTimecourse'].T,
                                                                                '''est['compsIdealAll']''',
                                                                                suite2p['compSpatial'])                                                  
        # extract the paired traces (simulated)                                                                 
        suite2p['simTraces'] = neur_act2[suite2p['pairs'][:,0],:]

        # extract the paired traces (estimated)
        suite2p['pairedTraces'] = suite2p['compTimecourse'][suite2p['pairs'][:,1],:]

    if calc_only[3]:
        # calculate pairs based on overlap and best temporal correlation
        est['corrvals'], est['pairs'], est['corrM'] = corrTimeAndSpace((neur_act2[sL]).T,
                                                                                est['estact'][0:-1,:].T,
                                                                                '''est['compsIdealAll']''',
                                                                       est['compsIdeal'])                                                  
        # extract the paired traces (simulated)                                                                 
        est['simTraces'] = neur_act2[est['pairs'][:,0],:]

        # extract the paired traces (estimated)
        est['pairedTraces'] = est['estact'][est['pairs'][:,1],:]
        
        
    return cnmf, pcaica, suite2p, est
