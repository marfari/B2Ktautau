from Configurables import DecayTreeTuple
from DecayTreeTuple.Configuration import *
from Configurables import LoKi__Hybrid__TupleTool
from Configurables import CheckPV
from Configurables import DaVinci
from PhysConf.Selections import FilterSelection
from PhysSelPython.Wrappers import AutomaticData

isMC = True
year = '2018'
pol = 'MagUp'
applySel = False

# redefine addBranches and setDescriptorTemplate because of weird bug in DV
def addBranches(self, branches):
    'Simplified adding of branches a little bit'
    'takes a dictionary of {branch: decay descriptor}, returns a dictionary of {branch: configurable instances}'
    from Configurables import TupleToolDecay
    if 'Branches' not in dir(self):
        raise TypeError('you\'re trying to add branches to something which doesn\'t support branching, ' + str(type(self)))
    if not isinstance(branches, type({})):
        raise TypeError('expected a dictionary of branches, got a ' + str(type(branches)) + ' instead')

    if self.Branches is None:
        self.Branches = {}

    instances = {}
    for branch in branches:
        # check for whitespace
        #for char in ""\t\r\s":
            #if char in branch:
                #raise NameError('You have tried to add a branch named \'' + branch + '\',which contains whitespace. This is not permitted.')
        self.Branches[branch] = branches[branch]
        self.addTool(TupleToolDecay, branch)
        instances[branch] = getattr(self, branch)

    return instances

def setDescriptorTemplate(self, template):
    'Bored of typing decay descriptors and adding carat symbols?'
    'Use some python string template magic to set your decay descriptor'
    'and define your branches all in one go without excess typing!'
    'from https://gitlab.cern.ch/snippets/147'
    if 'Decay' not in dir(self):
        raise TypeError('You\'re trying to set the decay descriptor of something that doesn\'t have one, ' + str(type(self)))
    if 'Branches' not in dir(self):
        raise TypeError('You\'re trying to define branches on something that doesn\'t support them, ' + str(type(self)))

    from string import Template
    # The argument 'template' is a Python string template
    # e.g. '${D}[D0 -> ${kaon}K- ${pion}pi+]CC'
    # Here ['D', 'kaon', 'pion'] are the branch names you want
    dd = Template(template)

    # This parses the temlate to get the list of branch names,
    # i.e. ['D', 'kaon', 'pion']
    particles = [y[1] if len(y[1]) else y[2] for y in dd.pattern.findall(dd.template) if len(y[1]) or len(y[2])]

    # To form the decay descriptor, we need to mark all the particles
    # except for the top-level particle, which is included by default
    mapping = {p: '^' for p in particles }
    mapping[particles[0]] = ''

    # Make the descriptor
    # '[D0 -> ^K- ^pi+]CC'
    self.Decay = dd.substitute(mapping)

    # Now make the branches
    branches = {}
    for p in particles:
        # Need a version of the descriptor where particle 'p' is marked but nothing else is.
        mapping = {q: '^' if p == q else '' for q in particles }
        branches[p] = dd.substitute(mapping)

    # Finally, add the branches to the DecayTreeTuple
    return self.addBranches(branches)

DecayTreeTuple.addBranches = addBranches
MCDecayTreeTuple.addBranches = addBranches
DecayTreeTuple.setDescriptorTemplate = setDescriptorTemplate
MCDecayTreeTuple.setDescriptorTemplate = setDescriptorTemplate

def addMInfo( branch ):
	MInfo = branch.addTupleTool("LoKi::Hybrid::TupleTool/M")
	MInfo.Variables = {
		"M12" : "M12",
		"M23" : "M23",
		"M13" : "M13"
	}

def addSubmassInfo( branch ):
	SubmassInfo = branch.addTupleTool("TupleToolSubMass")
	SubmassInfo.SetMax = 3

# Isolation information (stripping line output isolation variables)
def addIsoInfo(branch, stream, line, year):
	LoKi_Cone =  branch.Bp.addTupleTool("LoKi::Hybrid::TupleTool/KoKi_Cone")

	LoKi_Cone.Variables = {
		# Cone isolation variables (angle = 1.5)
		"CONEANGLE_15_B" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar1', 'CONEANGLE', -100)".format(stream,line),
		"CONEMULT_15_B" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar1', 'CONEMULT', -100)".format(stream,line),
		"CONEPTASYM_15_B" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar1', 'CONEPTASYM', -100)".format(stream,line),

		# Cone isolation variables (angle = 1.7)
		"CONEANGLE_17_B" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar2', 'CONEANGLE', -100)".format(stream,line),
		"CONEMULT_17_B" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar2', 'CONEMULT', -100)".format(stream,line),
		"CONEPTASYM_17_B" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar2', 'CONEPTASYM', -100)".format(stream,line),

		# Cone isolation variables (angle = 1.0)
		"CONEANGLE_10_B" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar3', 'CONEANGLE', -100)".format(stream,line),
		"CONEMULT_10_B" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar3', 'CONEMULT', -100)".format(stream,line),
		"CONEPTASYM_10_B" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar3', 'CONEPTASYM', -100)".format(stream,line),

	}

def addTopInfo( branch ):
    loki_top = branch.addTupleTool("LoKi::Hybrid::TupleTool/loki_top")
    loki_top.Preambulo += [
        "PFUNA = LoKi.Particles.PFunA"
    ]
    loki_top.Variables = {
        "PHI" : "PHI",
        "ETA" : "ETA",
        "MIPCHI2DV":"MIPCHI2DV(PRIMARY)",
        "BPVVDCHI2":"BPVVDCHI2",
        "BPVVD" : "BPVVD",
        "BPVDIRA":"BPVDIRA",
        "BPVLTIME":"BPVLTIME()",
        "BPVZ":"BPV(VZ)",
    	"VCHI2PDOF":"VFASPF(VCHI2PDOF)",
        # "AM" : "PFUNA(AM)"
    }

def addTrkInfo( branch ):
    loki_trk = branch.addTupleTool("LoKi::Hybrid::TupleTool/loki_trk")
    loki_trk.Variables = {
    	"PHI" : "PHI",
        "ETA" : "ETA",
        "MIPCHI2DV":"MIPCHI2DV(PRIMARY)",
        "TRCHI2DOF":"TRCHI2DOF",
        "TRGHOSTPROB":"TRGHOSTPROB"
    }

# DecayTreeFitter
def addKTTDTF(branch):
    # standard DTF
    dtf = branch.addTupleTool("TupleToolDecayTreeFitter/dtf")
    dtf.daughtersToConstrain = ['D0']
    dtf.UseFullTreeInName = True
    dtf.UpdateDaughters = True
    dtf.constrainToOriginVertex = True


def addTISTOS( branch ):
    ttb = branch.addTupleTool("TupleToolTISTOS")
    ttb.VerboseL0=True #if you want to fill info for each trigger line
    ttb.VerboseHlt1=True
    ttb.VerboseHlt2=True
    ttb.TriggerList=[
		# L0
		"L0HadronDecision",
		"L0MuonDecision",
		"L0ElectronDecision",
		"L0PhotonDecision",
		"L0CALODecision",
		# HLT1
		"Hlt1TrackMVADecision",
		"Hlt1TrackMVALooseDecision",
		"Hlt1TwoTrackMVADecision",
		"Hlt1TwoTrackMVALooseDecision",
		# "Hlt1CalibHighPTLowMultTrksDecision",
		# "Hlt1IncPhiDecision",
		# HLT2
		"Hlt2ForwardDecision",
		# "Hlt2XcMuXForTauB2XcFakeMuDecision",
		"Hlt2Topo2BodyDecision",
		"Hlt2Topo3BodyDecision",
		"Hlt2Topo4BodyDecision"
    ]
	
def addTupleToolL0Calo( branch ):
	l0_hcal_calo = branch.addTupleTool("TupleToolL0Calo")
	l0_hcal_calo.WhichCalo = "HCAL"

def addTools( DecayTreeTuple, stream, line, year ):
    #run modified DTF
    addKTTDTF(DecayTreeTuple.Bp) 
    #add some helpful variables for B0
    addTopInfo(DecayTreeTuple.Bp)
    #add isolation info
    addIsoInfo(DecayTreeTuple, stream, line, year)
    #add submasses
    addSubmassInfo(DecayTreeTuple.Bp)
    #add helpful info for final state tracks
    addTrkInfo(DecayTreeTuple.D0bar_K)
    addTrkInfo(DecayTreeTuple.D0bar_pi1)
    addTrkInfo(DecayTreeTuple.D0bar_pi2)
    addTrkInfo(DecayTreeTuple.D0bar_pi3)
    addTrkInfo(DecayTreeTuple.D0_K)
    addTrkInfo(DecayTreeTuple.D0_pi1)
    addTrkInfo(DecayTreeTuple.D0_pi2)
    addTrkInfo(DecayTreeTuple.D0_pi3)
    addTrkInfo(DecayTreeTuple.Kp)
    #add M12, M23, M13 for tau's
    addMInfo(DecayTreeTuple.D0bar)
    addMInfo(DecayTreeTuple.D0)
    # add trigger information
    addTISTOS(DecayTreeTuple.Bp)
    if isMC:
        addTupleToolL0Calo(DecayTreeTuple.D0bar_K)
        addTupleToolL0Calo(DecayTreeTuple.D0bar_pi1)
        addTupleToolL0Calo(DecayTreeTuple.D0bar_pi2)
        addTupleToolL0Calo(DecayTreeTuple.D0bar_pi3)
        addTupleToolL0Calo(DecayTreeTuple.D0_K)
        addTupleToolL0Calo(DecayTreeTuple.D0_pi1)
        addTupleToolL0Calo(DecayTreeTuple.D0_pi2)
        addTupleToolL0Calo(DecayTreeTuple.D0_pi3)
        addTupleToolL0Calo(DecayTreeTuple.Kp)

# Configure DaVinci
if isMC:
	# DaVinci().TupleFile = 'DVntuple_MC_'+year+'_'+pol+'.root'
	DaVinci().TupleFile = '/panfs/felician/BuD0D0K/DVntuple_MC_'+year+'_'+pol+'.root'
else:
	# DaVinci().TupleFile = 'DVntuple_data_'+year+'_'+pol+'.root'
	DaVinci().TupleFile = '/panfs/felician/BuD0D0K/DVntuple_data_'+year+'_'+pol+'.root'

# Stream and stripping line we want to use
if isMC:
	stream = 'AllStreams' 
else:	
	stream = 'Bhadron'
line = 'B2D0D0KD02K3PiD02K3PiBeauty2CharmLine'

if isMC:
	DaVinci().Simulation = True

	# CondDBtag / DDDBtag specify the exact detector conditions with which the MC was generated. They are specified in the downloaded DST file.
	# if year == '2016':
	# 	DaVinci().DDDBtag = "dddb-20220927-2016" 
	# 	if pol == 'MagUp':
	# 		DaVinci().CondDBtag = "sim-20201113-6-vc-mu100-Sim10"  
	# 	elif pol == 'MagDown':
	# 		DaVinci().CondDBtag = "sim-20201113-6-vc-md100-Sim10"
	if year == '2017':
		DaVinci().DDDBtag = "dddb-20170721-3"
		if pol == 'MagUp':
			DaVinci().CondDBtag = "sim-20190430-1-vc-mu100"
		elif pol == 'MagDown':
			DaVinci().CondDBtag = "sim-20190430-1-vc-md100"
	elif year == '2018':
		DaVinci().DDDBtag = "dddb-20170721-3"
		if pol == 'MagUp':
			DaVinci().CondDBtag = "sim-20190430-vc-mu100" 
		elif pol == 'MagDown':
			DaVinci().CondDBtag = "sim-20190430-vc-md100"
else:
	DaVinci().Simulation = False
DaVinci().InputType = 'MDST'
DaVinci().RootInTES = '/Event/{0}'.format(stream) # for MDST only

DaVinci().PrintFreq = 1000
DaVinci().DataType = year 

DaVinci().Lumi = not DaVinci().Simulation # Only ask for luminosity information when not using simulated data
DaVinci().EvtMax = -1 # how many events to run over (-1 = runs over all events)

# Add smearing to MC and momentum scaling to data
if isMC:
	# Apply momentum smearing
	from Configurables import TrackSmearState as SMEAR
	smear = SMEAR('StateSmear')
	DaVinci().UserAlgorithms += [smear]
else:
	# Apply momentum scaling to data
	from Configurables import TrackScaleState as SCALE
	from Configurables import CondDB
	CondDB(LatestGlobalTagByDataType=year) # sets the DDDB and CondDB for data to the latest recommended tags for that year
	scaler = SCALE('StateScale')
	DaVinci().UserAlgorithms += [scaler]

# Create an ntuple to capture B+ decays from the StrippingLine line
dtt = DecayTreeTuple('ntuple')
mcdtt = MCDecayTreeTuple('mc_ntuple')

# Loose preselections to reduce the size of the tuples (in data)
# These (rectangular) selections are applied to MC offline
# filter code (will be updated)
filtercode = ""

if applySel:
	if isMC:
		b_strip = AutomaticData('/Event/{0}/Phys/{1}/Particles'.format(stream, line))
	if not isMC:
		b_strip = AutomaticData('/Phys/{0}/Particles'.format(line)) 

	b_Filter    = FilterSelection("b_Filter", b_strip, Code = filtercode)
	dtt.Inputs = [b_Filter.outputLocation()]
else:
	dtt.Inputs = ['Phys/{0}/Particles'.format(line)]

dtt.setDescriptorTemplate('${Bp}[B+ -> ${Kp}K+ ${D0bar}(D0 -> ${D0bar_K}K+ ${D0bar_pi1}pi+ ${D0bar_pi2}pi- ${D0bar_pi3}pi-) ${D0}(D0 -> ${D0_K}K+ ${D0_pi1}pi+ ${D0_pi2}pi- ${D0_pi3}pi-)]CC')

# Generated (Acceptance cut: DaughtersInLHCb)
mcdtt.Decay = '[B+ => ^(D~0 => ^K+ ^pi+ ^pi- ^pi-)  ^(D0 => ^K- ^pi- ^pi+ ^pi+) K+]CC' 

# TUPLE TOOLS
ToolList = ['TupleToolTrackInfo', 'TupleToolVtxIsoln', 'TupleToolPrimaries', 'TupleToolCorrectedMass', 'TupleToolPropertime', 'TupleToolAngles', 'TupleToolL0Calo', 'TupleToolRecoStats']
dtt.ToolList += ToolList
if isMC:
	dtt.addTupleTool('TupleToolMCTruth')
	dtt.addTupleTool('TupleToolMCBackgroundInfo')

# Adds several tupletools, including DTF
addTools(dtt, stream, line, year)

GENToolList = ['MCTupleToolKinematic', 'MCTupleToolHierarchy', 'TupleToolEventInfo']
if isMC:
	mcdtt.ToolList += GENToolList

# Stripping filter 
from PhysConf.Filters import LoKi_Filters
fltrs = LoKi_Filters (
    STRIP_Code = "HLT_PASS_RE('StrippingB2D0D0KD02K3PiD02K3PiBeauty2CharmLineDecision')"
)
if not isMC: # don't apply to MC because it will affect the MCDecayTree
	DaVinci().EventPreFilters = fltrs.filters('Filters')

# PV checker -> TupleToolPrimaries requires there to be a PV in the event, otherwise it throws an error
DaVinci().UserAlgorithms += [CheckPV()]

if applySel:
	DaVinci().UserAlgorithms += [b_Filter]

DaVinci().UserAlgorithms += [dtt]
if isMC:
	DaVinci().UserAlgorithms += [mcdtt]

# Use the local input data
from GaudiConf import IOHelper
IOHelper('ROOT').inputFiles([
    '/panfs/felician/BuD0D0K/DST_MC_2018_MagUp/00116385_00000002_7.AllStreams.mdst',
    '/panfs/felician/BuD0D0K/DST_MC_2018_MagUp/00116385_00000003_7.AllStreams.mdst',
    '/panfs/felician/BuD0D0K/DST_MC_2018_MagUp/00116385_00000012_7.AllStreams.mdst',
    '/panfs/felician/BuD0D0K/DST_MC_2018_MagUp/00116385_00000007_7.AllStreams.mdst'
	# '/panfs/felician/BuD0D0K/2018/DST_data_MagUp/00076476_00010124_1.bhadron.mdst'
], clear=True)