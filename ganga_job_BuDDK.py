isMC = False
year = '2016'
pol = 'MagDown'

j = Job(name='D6_Down_DDK_noSel')
myApp = GaudiExec()
myApp.directory = "DaVinciDev"
j.application = myApp
j.application.options = ['ntuple_options_BuDDK.py']
j.application.platform = 'x86_64-centos7-gcc9-opt'

if isMC:
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

j.inputdata = data[:10]
j.backend = Dirac()
j.splitter = SplitByFiles(filesPerJob=5) # ~20 for MC; 5 for data

j.outputfiles = [DiracFile('*.root'),LocalFile('stdout')]

j.submit()
