isMC = False
year = '2016'
pol = 'MagDown'

j = Job(name='2016 MagDown Data')
myApp = GaudiExec()
myApp.directory = "DaVinciDev"
j.application = myApp
j.application.options = ['ntuple_options.py']
j.application.platform = 'x86_64-centos7-gcc9-opt'

if not isMC:
    if((year == '2011') and (pol == 'MagDown')):
        bkPath = '/LHCb/Collision11/Beam3500GeV-VeloClosed-MagDown/Real Data/Reco14/Stripping21r1p2/90000000/BHADRON.MDST'
    elif((year == '2011') and (pol == 'MagUp')):
        bkPath = '/LHCb/Collision11/Beam3500GeV-VeloClosed-MagUp/Real Data/Reco14/Stripping21r1p2/90000000/BHADRON.MDST'
    elif((year == '2012') and (pol == 'MagDown')):
        bkPath = '/LHCb/Collision12/Beam4000GeV-VeloClosed-MagDown/Real Data/Reco14/Stripping21r0p2/90000000/BHADRON.MDST'
    elif((year == '2012') and (pol == 'MagUp')):
        bkPath = '/LHCb/Collision12/Beam4000GeV-VeloClosed-MagUp/Real Data/Reco14/Stripping21r0p2/90000000/BHADRON.MDST'
    elif((year == '2015') and (pol == "MagDown")):
        bkPath = '/LHCb/Collision15/Beam6500GeV-VeloClosed-MagDown/Real Data/Reco15a/Stripping24r2/90000000/BHADRON.MDST'  
    elif((year == '2015') and (pol == 'MagUp')):
        bkPath = '/LHCb/Collision15/Beam6500GeV-VeloClosed-MagUp/Real Data/Reco15a/Stripping24r2/90000000/BHADRON.MDST'
    elif((year == '2016') and (pol == 'MagDown')):
         bkPath = '/LHCb/Collision16/Beam6500GeV-VeloClosed-MagDown/Real Data/Reco16/Stripping28r2/90000000/BHADRON.MDST' 
    elif((year == '2016') and (pol == 'MagUp')):
        bkPath = '/LHCb/Collision16/Beam6500GeV-VeloClosed-MagUp/Real Data/Reco16/Stripping28r2/90000000/BHADRON.MDST' 
    elif((year == '2017') and (pol == 'MagDown')):
        bkPath = '/LHCb/Collision17/Beam6500GeV-VeloClosed-MagDown/Real Data/Reco17/Stripping29r2p1/90000000/BHADRON.MDST'
    elif((year == '2017' and (pol == 'MagUp'))):
        bkPath = '/LHCb/Collision17/Beam6500GeV-VeloClosed-MagUp/Real Data/Reco17/Stripping29r2p1/90000000/BHADRON.MDST' 
    elif((year == '2018') and (pol == 'MagDown')):
        bkPath = '/LHCb/Collision18/Beam6500GeV-VeloClosed-MagDown/Real Data/Reco18/Stripping34r0p1/90000000/BHADRON.MDST' 
    elif((year == '2018') and (pol == 'MagUp')):
        bkPath = '/LHCb/Collision18/Beam6500GeV-VeloClosed-MagUp/Real Data/Reco18/Stripping34r0p1/90000000/BHADRON.MDST' 
data  = BKQuery(bkPath, dqflag=['OK']).getDataset()

j.inputdata = data
j.backend = Dirac()
j.splitter = SplitByFiles(filesPerJob=20)

if isMC:
    j.outputfiles = [LocalFile('DVntuple_MC_'+year+'_'+pol+'.root')] # change, put in grid DiracFile
else:
    #j.outputfiles = [LocalFile('DVntuple_data_'+year+'_'+pol+'.root')]
    j.outputfiles = [DiracFile('*.root'),LocalFile('stdout')]

j.submit()

