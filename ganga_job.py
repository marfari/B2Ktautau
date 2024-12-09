isMC = True
year = '2018'
pol = 'MagDown'
name = 'BuDD'

if(isMC):
    typeName = 'MC'
else:
    typeName = 'Data'

j = Job(name=typeName+'_'+year+'_'+pol+'_'+name)
myApp = GaudiExec()
myApp.directory = "DaVinciDev_v46r10"
j.application = myApp

if(name == 'BuKtautau'):
    j.application.options = ['ntuple_options.py']
elif(name == 'BuDDKp'):
    j.application.options = ['ntuple_options_BuDDKp_cocktail.py']
elif(name == 'BdDDKp'):
    j.application.options = ['ntuple_options_BdDDKp_cocktail.py']
elif(name == 'BsDDKp'):
    j.application.options = ['ntuple_options_BsDDKp_cocktail.py']
elif(name == 'BuDDK0'):
    j.application.options = ['ntuple_options_BuDDK0_cocktail.py']
elif(name == 'BdDDK0'):
    j.application.options = ['ntuple_options_BdDDK0_cocktail.py']
elif(name == 'BuDD'):
    j.application.options = ['ntuple_options_BuDD_cocktail.py']
elif(name == 'BdDD'):
    j.application.options = ['ntuple_options_BdDD_cocktail.py']
elif(name == 'BsDD'):
    j.application.options = ['ntuple_options_BsDD_cocktail.py']

j.application.platform = 'x86_64_v2-el9-gcc12-opt' 
# j.inputfiles += [LocalFile("/panfs/felician/MLP_weights/KTauTau_MLP_Train_taup_PX/dataset/weights/TMVARegression_taup_TRUEP_X_MLP.weights.xml"),
#                 LocalFile("/panfs/felician/MLP_weights/KTauTau_MLP_Train_taup_PX/dataset/weights/TMVARegression_taup_TRUEP_X_MLP_fold1.weights.xml"),
#                 LocalFile("/panfs/felician/MLP_weights/KTauTau_MLP_Train_taup_PX/dataset/weights/TMVARegression_taup_TRUEP_X_MLP_fold2.weights.xml"),
#                 LocalFile("/panfs/felician/MLP_weights/KTauTau_MLP_Train_taup_PY/dataset/weights/TMVARegression_taup_TRUEP_Y_MLP.weights.xml"),
#                 LocalFile("/panfs/felician/MLP_weights/KTauTau_MLP_Train_taup_PY/dataset/weights/TMVARegression_taup_TRUEP_Y_MLP_fold1.weights.xml"),
#                 LocalFile("/panfs/felician/MLP_weights/KTauTau_MLP_Train_taup_PY/dataset/weights/TMVARegression_taup_TRUEP_Y_MLP_fold2.weights.xml"),
#                 LocalFile("/panfs/felician/MLP_weights/KTauTau_MLP_Train_taup_PZ/dataset/weights/TMVARegression_taup_TRUEP_Z_MLP.weights.xml"),
#                 LocalFile("/panfs/felician/MLP_weights/KTauTau_MLP_Train_taup_PZ/dataset/weights/TMVARegression_taup_TRUEP_Z_MLP_fold1.weights.xml"),
#                 LocalFile("/panfs/felician/MLP_weights/KTauTau_MLP_Train_taup_PZ/dataset/weights/TMVARegression_taup_TRUEP_Z_MLP_fold2.weights.xml"),
#                 LocalFile("/panfs/felician/MLP_weights/KTauTau_MLP_Train_taum_PX/dataset/weights/TMVARegression_taum_TRUEP_X_MLP.weights.xml"),
#                 LocalFile("/panfs/felician/MLP_weights/KTauTau_MLP_Train_taum_PX/dataset/weights/TMVARegression_taum_TRUEP_X_MLP_fold1.weights.xml"),
#                 LocalFile("/panfs/felician/MLP_weights/KTauTau_MLP_Train_taum_PX/dataset/weights/TMVARegression_taum_TRUEP_X_MLP_fold2.weights.xml"),
#                 LocalFile("/panfs/felician/MLP_weights/KTauTau_MLP_Train_taum_PY/dataset/weights/TMVARegression_taum_TRUEP_Y_MLP.weights.xml"),
#                 LocalFile("/panfs/felician/MLP_weights/KTauTau_MLP_Train_taum_PY/dataset/weights/TMVARegression_taum_TRUEP_Y_MLP_fold1.weights.xml"),
#                 LocalFile("/panfs/felician/MLP_weights/KTauTau_MLP_Train_taum_PY/dataset/weights/TMVARegression_taum_TRUEP_Y_MLP_fold2.weights.xml"),
#                 LocalFile("/panfs/felician/MLP_weights/KTauTau_MLP_Train_taum_PZ/dataset/weights/TMVARegression_taum_TRUEP_Z_MLP.weights.xml"),
#                 LocalFile("/panfs/felician/MLP_weights/KTauTau_MLP_Train_taum_PZ/dataset/weights/TMVARegression_taum_TRUEP_Z_MLP_fold1.weights.xml"),
#                 LocalFile("/panfs/felician/MLP_weights/KTauTau_MLP_Train_taum_PZ/dataset/weights/TMVARegression_taum_TRUEP_Z_MLP_fold2.weights.xml") 
#                 ]

if isMC:
    if(name == 'BuKtautau'):
        if((year == '2016') and (pol == 'MagDown')):
            bkPath = '/MC/2016/Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8/Sim10c/Trig0x6139160F/Reco16/Turbo03a/Stripping28r2p2Filtered/12201011/B2KTAUTAU.STRIP.DST'
        elif((year == '2016') and (pol == 'MagUp')):
            bkPath = '/MC/2016/Beam6500GeV-2016-MagUp-Nu1.6-25ns-Pythia8/Sim10c/Trig0x6139160F/Reco16/Turbo03a/Stripping28r2p2Filtered/12201011/B2KTAUTAU.STRIP.DST'
        elif((year == '2017') and (pol == 'MagDown')):
            bkPath = '/MC/2017/Beam6500GeV-2017-MagDown-Nu1.6-25ns-Pythia8/Sim10c/Trig0x62661709/Reco17/Turbo04a-WithTurcal/Stripping29r2p3Filtered/12201011/B2KTAUTAU.STRIP.DST'
        elif((year == '2017' and (pol == 'MagUp'))):
            bkPath = '/MC/2017/Beam6500GeV-2017-MagUp-Nu1.6-25ns-Pythia8/Sim10c/Trig0x62661709/Reco17/Turbo04a-WithTurcal/Stripping29r2p3Filtered/12201011/B2KTAUTAU.STRIP.DST'
        elif((year == '2018') and (pol == 'MagDown')):
            bkPath = '/MC/2018/Beam6500GeV-2018-MagDown-Nu1.6-25ns-Pythia8/Sim10c/Trig0x617d18a4/Reco18/Turbo05-WithTurcal/Stripping34r0p3Filtered/12201011/B2KTAUTAU.STRIP.DST'
        elif((year == '2018') and (pol == 'MagUp')):
            bkPath = '/MC/2018/Beam6500GeV-2018-MagUp-Nu1.6-25ns-Pythia8/Sim10c/Trig0x617d18a4/Reco18/Turbo05-WithTurcal/Stripping34r0p3Filtered/12201011/B2KTAUTAU.STRIP.DST'
    elif(name == 'BuDDKp'):
        if((year == '2016') and (pol == 'MagDown')):
            bkPath = '/MC/2016/Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8/Sim10d/Trig0x6139160F/Reco16/Turbo03a/Stripping28r2p2Filtered/12693500/B2KTAUTAU.STRIP.DST'
        elif((year == '2016') and (pol == 'MagUp')):
            bkPath = '/MC/2016/Beam6500GeV-2016-MagUp-Nu1.6-25ns-Pythia8/Sim10d/Trig0x6139160F/Reco16/Turbo03a/Stripping28r2p2Filtered/12693500/B2KTAUTAU.STRIP.DST'
        elif((year == '2017') and (pol == 'MagDown')):
            bkPath = '/MC/2017/Beam6500GeV-2017-MagDown-Nu1.6-25ns-Pythia8/Sim10d/Trig0x62661709/Reco17/Turbo04a-WithTurcal/Stripping29r2p3Filtered/12693500/B2KTAUTAU.STRIP.DST'
        elif((year == '2017' and (pol == 'MagUp'))):
            bkPath = '/MC/2017/Beam6500GeV-2017-MagUp-Nu1.6-25ns-Pythia8/Sim10d/Trig0x62661709/Reco17/Turbo04a-WithTurcal/Stripping29r2p3Filtered/12693500/B2KTAUTAU.STRIP.DST'
        elif((year == '2018') and (pol == 'MagDown')):
            bkPath = '/MC/2018/Beam6500GeV-2018-MagDown-Nu1.6-25ns-Pythia8/Sim10d/Trig0x617d18a4/Reco18/Turbo05-WithTurcal/Stripping34r0p3Filtered/12693500/B2KTAUTAU.STRIP.DST'
        elif((year == '2018') and (pol == 'MagUp')):
            bkPath = '/MC/2018/Beam6500GeV-2018-MagUp-Nu1.6-25ns-Pythia8/Sim10d/Trig0x617d18a4/Reco18/Turbo05-WithTurcal/Stripping34r0p3Filtered/12693500/B2KTAUTAU.STRIP.DST'
    elif(name == 'BdDDKp'):
        if((year == '2016') and (pol == 'MagDown')):
            bkPath = '/MC/2016/Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8/Sim10d/Trig0x6139160F/Reco16/Turbo03a/Stripping28r2p2Filtered/11294500/B2KTAUTAU.STRIP.DST'
        elif((year == '2016') and (pol == 'MagUp')):
            bkPath = '/MC/2016/Beam6500GeV-2016-MagUp-Nu1.6-25ns-Pythia8/Sim10d/Trig0x6139160F/Reco16/Turbo03a/Stripping28r2p2Filtered/11294500/B2KTAUTAU.STRIP.DST'
        elif((year == '2017') and (pol == 'MagDown')):
            bkPath = '/MC/2017/Beam6500GeV-2017-MagDown-Nu1.6-25ns-Pythia8/Sim10d/Trig0x62661709/Reco17/Turbo04a-WithTurcal/Stripping29r2p3Filtered/11294500/B2KTAUTAU.STRIP.DST'
        elif((year == '2017' and (pol == 'MagUp'))):
            bkPath = '/MC/2017/Beam6500GeV-2017-MagUp-Nu1.6-25ns-Pythia8/Sim10d/Trig0x62661709/Reco17/Turbo04a-WithTurcal/Stripping29r2p3Filtered/11294500/B2KTAUTAU.STRIP.DST'
        elif((year == '2018') and (pol == 'MagDown')):
            bkPath = '/MC/2018/Beam6500GeV-2018-MagDown-Nu1.6-25ns-Pythia8/Sim10d/Trig0x617d18a4/Reco18/Turbo05-WithTurcal/Stripping34r0p3Filtered/11294500/B2KTAUTAU.STRIP.DST'
        elif((year == '2018') and (pol == 'MagUp')):
            bkPath = '/MC/2018/Beam6500GeV-2018-MagUp-Nu1.6-25ns-Pythia8/Sim10d/Trig0x617d18a4/Reco18/Turbo05-WithTurcal/Stripping34r0p3Filtered/11294500/B2KTAUTAU.STRIP.DST'
    elif(name == 'BsDDKp'):
        if((year == '2016') and (pol == 'MagDown')):
            bkPath = '/MC/2016/Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8/Sim10d/Trig0x6139160F/Reco16/Turbo03a/Stripping28r2p2Filtered/13297500/B2KTAUTAU.STRIP.DST'
        elif((year == '2016') and (pol == 'MagUp')):
            bkPath = '/MC/2016/Beam6500GeV-2016-MagUp-Nu1.6-25ns-Pythia8/Sim10d/Trig0x6139160F/Reco16/Turbo03a/Stripping28r2p2Filtered/13297500/B2KTAUTAU.STRIP.DST'
        elif((year == '2017') and (pol == 'MagDown')):
            bkPath = '/MC/2017/Beam6500GeV-2017-MagDown-Nu1.6-25ns-Pythia8/Sim10d/Trig0x62661709/Reco17/Turbo04a-WithTurcal/Stripping29r2p3Filtered/13297500/B2KTAUTAU.STRIP.DST'
        elif((year == '2017' and (pol == 'MagUp'))):
            bkPath = '/MC/2017/Beam6500GeV-2017-MagUp-Nu1.6-25ns-Pythia8/Sim10d/Trig0x62661709/Reco17/Turbo04a-WithTurcal/Stripping29r2p3Filtered/13297500/B2KTAUTAU.STRIP.DST'
        elif((year == '2018') and (pol == 'MagDown')):
            bkPath = '/MC/2018/Beam6500GeV-2018-MagDown-Nu1.6-25ns-Pythia8/Sim10d/Trig0x617d18a4/Reco18/Turbo05-WithTurcal/Stripping34r0p3Filtered/13297500/B2KTAUTAU.STRIP.DST'
        elif((year == '2018') and (pol == 'MagUp')):
            bkPath = '/MC/2018/Beam6500GeV-2018-MagUp-Nu1.6-25ns-Pythia8/Sim10d/Trig0x617d18a4/Reco18/Turbo05-WithTurcal/Stripping34r0p3Filtered/13297500/B2KTAUTAU.STRIP.DST'
    elif(name == 'BuDDK0'):
        if((year == '2016') and (pol == 'MagDown')):
            bkPath = '/MC/2016/Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8/Sim10d/Trig0x6139160F/Reco16/Turbo03a/Stripping28r2p2Filtered/12293500/B2KTAUTAU.STRIP.DST'
        elif((year == '2016') and (pol == 'MagUp')):
            bkPath = '/MC/2016/Beam6500GeV-2016-MagUp-Nu1.6-25ns-Pythia8/Sim10d/Trig0x6139160F/Reco16/Turbo03a/Stripping28r2p2Filtered/12293500/B2KTAUTAU.STRIP.DST'
        elif((year == '2017') and (pol == 'MagDown')):
            bkPath = '/MC/2017/Beam6500GeV-2017-MagDown-Nu1.6-25ns-Pythia8/Sim10d/Trig0x62661709/Reco17/Turbo04a-WithTurcal/Stripping29r2p3Filtered/12293500/B2KTAUTAU.STRIP.DST'
        elif((year == '2017' and (pol == 'MagUp'))):
            bkPath = '/MC/2017/Beam6500GeV-2017-MagUp-Nu1.6-25ns-Pythia8/Sim10d/Trig0x62661709/Reco17/Turbo04a-WithTurcal/Stripping29r2p3Filtered/12293500/B2KTAUTAU.STRIP.DST'
        elif((year == '2018') and (pol == 'MagDown')):
            bkPath = '/MC/2018/Beam6500GeV-2018-MagDown-Nu1.6-25ns-Pythia8/Sim10d/Trig0x617d18a4/Reco18/Turbo05-WithTurcal/Stripping34r0p3Filtered/12293500/B2KTAUTAU.STRIP.DST'
        elif((year == '2018') and (pol == 'MagUp')):
            bkPath = '/MC/2018/Beam6500GeV-2018-MagUp-Nu1.6-25ns-Pythia8/Sim10d/Trig0x617d18a4/Reco18/Turbo05-WithTurcal/Stripping34r0p3Filtered/12293500/B2KTAUTAU.STRIP.DST'
    elif(name == 'BuDD'):
        if((year == '2016') and (pol == 'MagDown')):
            bkPath = '/MC/2016/Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8/Sim10d/Trig0x6139160F/Reco16/Turbo03a/Stripping28r2p2Filtered/12696700/B2KTAUTAU.STRIP.DST'
        elif((year == '2016') and (pol == 'MagUp')):
            bkPath = '/MC/2016/Beam6500GeV-2016-MagUp-Nu1.6-25ns-Pythia8/Sim10d/Trig0x6139160F/Reco16/Turbo03a/Stripping28r2p2Filtered/12696700/B2KTAUTAU.STRIP.DST'
        elif((year == '2017') and (pol == 'MagDown')):
            bkPath = '/MC/2017/Beam6500GeV-2017-MagDown-Nu1.6-25ns-Pythia8/Sim10d/Trig0x62661709/Reco17/Turbo04a-WithTurcal/Stripping29r2p3Filtered/12696700/B2KTAUTAU.STRIP.DST'
        elif((year == '2017' and (pol == 'MagUp'))):
            bkPath = '/MC/2017/Beam6500GeV-2017-MagUp-Nu1.6-25ns-Pythia8/Sim10d/Trig0x62661709/Reco17/Turbo04a-WithTurcal/Stripping29r2p3Filtered/12696700/B2KTAUTAU.STRIP.DST'
        elif((year == '2018') and (pol == 'MagDown')):
            bkPath = '/MC/2018/Beam6500GeV-2018-MagDown-Nu1.6-25ns-Pythia8/Sim10d/Trig0x617d18a4/Reco18/Turbo05-WithTurcal/Stripping34r0p3Filtered/12696700/B2KTAUTAU.STRIP.DST'
        elif((year == '2018') and (pol == 'MagUp')):
            bkPath = '/MC/2018/Beam6500GeV-2018-MagUp-Nu1.6-25ns-Pythia8/Sim10d/Trig0x617d18a4/Reco18/Turbo05-WithTurcal/Stripping34r0p3Filtered/12696700/B2KTAUTAU.STRIP.DST'
else:
    if((year == '2016') and (pol == 'MagDown')):
         bkPath = '/LHCb/Collision16/Beam6500GeV-VeloClosed-MagDown/Real Data/Reco16/Stripping28r2p2/90000000/BHADRON.MDST' 
    elif((year == '2016') and (pol == 'MagUp')):
        bkPath = '/LHCb/Collision16/Beam6500GeV-VeloClosed-MagUp/Real Data/Reco16/Stripping28r2p2/90000000/BHADRON.MDST' 
    elif((year == '2017') and (pol == 'MagDown')):
        bkPath = '/LHCb/Collision17/Beam6500GeV-VeloClosed-MagDown/Real Data/Reco17/Stripping29r2p3/90000000/BHADRON.MDST'
    elif((year == '2017' and (pol == 'MagUp'))):
        bkPath = '/LHCb/Collision17/Beam6500GeV-VeloClosed-MagUp/Real Data/Reco17/Stripping29r2p3/90000000/BHADRON.MDST' 
    elif((year == '2018') and (pol == 'MagDown')):
        bkPath = '/LHCb/Collision18/Beam6500GeV-VeloClosed-MagDown/Real Data/Reco18/Stripping34r0p3/90000000/BHADRON.MDST' 
    elif((year == '2018') and (pol == 'MagUp')):
        bkPath = '/LHCb/Collision18/Beam6500GeV-VeloClosed-MagUp/Real Data/Reco18/Stripping34r0p3/90000000/BHADRON.MDST' 
data  = BKQuery(bkPath, dqflag=['OK']).getDataset()

j.inputdata = data
j.backend = Dirac()
if isMC:
    j.splitter = SplitByFiles(filesPerJob=20) # 20
else:
    j.splitter = SplitByFiles(filesPerJob=10) 

j.outputfiles = [DiracFile('*.root'),LocalFile('stdout')]

j.submit()

