isMC = False
year = '2016'
pol = 'MagUp'

j = Job(name='B2DDK 2016 MagUp')
#myApp = prepareGaudiExec('DaVinci','v46r4',myPath='.')
myApp = GaudiExec()
myApp.directory = "./DaVinciDev_v46r4"
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
data  = BKQuery(bkPath, dqflag=['OK']).getDataset()

j.inputdata = data[0:20]
j.backend = Dirac()
j.splitter = SplitByFiles(filesPerJob=5)

if isMC:
    j.outputfiles = [LocalFile('DVntuple_MC_'+year+'_'+pol+'.root')]
else:
    j.outputfiles = [LocalFile('DVntuple_data_'+year+'_'+pol+'.root')]

j.submit()
