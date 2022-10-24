isMC = False
year = '2016'
pol = 'MagUp'

j = Job(name='B2DDK 2015 MagUp')
#myApp = prepareGaudiExec('DaVinci','v46r4',myPath='.')
myApp = GaudiExec()
myApp.directory = "DaVinciDev_v46r4"
j.application = myApp
j.application.options = ['ntuple_options_BDDK.py']
j.application.platform = 'x86_64_v2-centos7-gcc11-opt'

if not isMC:
    if((year == '2011') and (pol == 'MagDown')):
        bkPath = '/LHCb/Collision11/Beam3500GeV-VeloClosed-MagDown/Real Data/Reco14/Stripping21r1/90000000/BHADRON.MDST'
    elif((year == '2011') and (pol == 'MagUp')):
        bkPath = '/LHCb/Collision11/Beam3500GeV-VeloClosed-MagUp/Real Data/Reco14/Stripping21r1/90000000/BHADRON.MDST'
    elif((year == '2012') and (pol == 'MagDown')):
        bkPath = '/LHCb/Collision12/Beam4000GeV-VeloClosed-MagDown/Real Data/Reco14/Stripping21/90000000/BHADRON.MDST'
    elif((year == '2012') and (pol == 'MagUp')):
        bkPath = '/LHCb/Collision12/Beam4000GeV-VeloClosed-MagUp/Real Data/Reco14/Stripping21/90000000/BHADRON.MDST'
    elif((year == '2015') and (pol == "MagDown")):
        bkPath = '/LHCb/Collision15/Beam6500GeV-VeloClosed-MagDown/Real Data/Reco15a/Stripping24r1/90000000/BHADRON.MDST'  
    elif((year == '2015') and (pol == 'MagUp')):
        bkPath = '/LHCb/Collision15/Beam6500GeV-VeloClosed-MagUp/Real Data/Reco15a/Stripping24r1/90000000/BHADRON.MDST'
    elif((year == '2016') and (pol == 'MagDown')):
         bkPath = '/LHCb/Collision16/Beam6500GeV-VeloClosed-MagDown/Real Data/Reco16/Stripping28r1/90000000/BHADRON.MDST' 
    elif((year == '2016') and (pol == 'MagUp')):
        bkPath = '/LHCb/Collision16/Beam6500GeV-VeloClosed-MagUp/Real Data/Reco16/Stripping28r1/90000000/BHADRON.MDST' 
    elif((year == '2017') and (pol == 'MagDown')):
        bkPath = '/LHCb/Collision17/Beam6500GeV-VeloClosed-MagDown/Real Data/Reco17/Stripping29r2/90000000/BHADRON.MDST'
    elif((year == '2017' and (pol == 'MagUp'))):
        bkPath = '/LHCb/Collision17/Beam6500GeV-VeloClosed-MagUp/Real Data/Reco17/Stripping29r2/90000000/BHADRON.MDST' 
    elif((year == '2018') and (pol == 'MagDown')):
        bkPath = '/LHCb/Collision18/Beam6500GeV-VeloClosed-MagDown/Real Data/Reco18/Stripping34/90000000/BHADRON.MDST' 
    elif((year == '2018') and (pol == 'MagUp')):
        bkPath = '/LHCb/Collision18/Beam6500GeV-VeloClosed-MagUp/Real Data/Reco18/Stripping34/90000000/BHADRON.MDST' 
elif isMC:
    if((year == '2015') and (pol == 'MagDown')):
        bkPath = '/MC/2015/Beam6500GeV-2015-MagDown-Nu1.6-25ns-Pythia8/Sim09c/Trig0x411400a2/Reco15a/Turbo02/Stripping24r1NoPrescalingFlagged/12197066/ALLSTREAMS.MDST'
    elif((year == '2015') and (pol == 'MagUp')):
        bkPath = '/MC/2015/Beam6500GeV-2015-MagUp-Nu1.6-25ns-Pythia8/Sim09c/Trig0x411400a2/Reco15a/Turbo02/Stripping24r1NoPrescalingFlagged/12197046/ALLSTREAMS.MDST'
    elif((year == '2016') and (pol == 'MagDown')):
        bkPath = '/MC/2016/Beam6500GeV-2016-MagDown-Nu1.6-25ns-Pythia8/Sim09c/Trig0x6138160F/Reco16/Turbo03/Stripping28r1NoPrescalingFlagged/12197046/ALLSTREAMS.MDST'
    elif((year == '2016') and (pol == 'MagUp')):
        bkPath = '/MC/2016/Beam6500GeV-2016-MagUp-Nu1.6-25ns-Pythia8/Sim09c/Trig0x6138160F/Reco16/Turbo03/Stripping28r1NoPrescalingFlagged/12197046/ALLSTREAMS.MDST'
    elif((year == '2017') and (pol == 'MagDown')):
        bkPath = '/MC/2017/Beam6500GeV-2017-MagDown-Nu1.6-25ns-Pythia8/Sim09e-ReDecay01/Trig0x62661709/Reco17/Turbo04a-WithTurcal/Stripping29r2NoPrescalingFlagged/12197046/ALLSTREAMS.MDST'
    elif((year == '2017') and (pol == 'MagUp')):
        bkPath = '/MC/2017/Beam6500GeV-2017-MagUp-Nu1.6-25ns-Pythia8/Sim09e-ReDecay01/Trig0x62661709/Reco17/Turbo04a-WithTurcal/Stripping29r2NoPrescalingFlagged/12197046/ALLSTREAMS.MDST'
    elif((year == '2018') and (pol == 'MagDown')):
        blPath = '/MC/2018/Beam6500GeV-2018-MagDown-Nu1.6-25ns-Pythia8/Sim09f-ReDecay01/Trig0x617d18a4/Reco18/Turbo05-WithTurcal/Stripping34NoPrescalingFlagged/12197046/ALLSTREAMS.MDST'
    elif((year == '2018') and (pol == 'MagUp')):
        bkPath = '/MC/2018/Beam6500GeV-2018-MagUp-Nu1.6-25ns-Pythia8/Sim09f-ReDecay01/Trig0x617d18a4/Reco18/Turbo05-WithTurcal/Stripping34NoPrescalingFlagged/12197046/ALLSTREAMS.MDST'
data  = BKQuery(bkPath, dqflag=['OK']).getDataset()

j.inputdata = data[0:1]
j.backend = Dirac()
j.splitter = SplitByFiles(filesPerJob=5)

if isMC:
    j.outputfiles = [LocalFile('DVntuple_MC_'+year+'_'+pol+'.root')]
else:
    j.outputfiles = [LocalFile('DVntuple_data_'+year+'_'+pol+'.root')]

j.submit()
