from Configurables import DecayTreeTuple
from DecayTreeTuple.Configuration import *
from Configurables import LoKi__Hybrid__TupleTool
from Configurables import CheckPV
from Configurables import DaVinci
from PhysConf.Selections import FilterSelection
from PhysSelPython.Wrappers import AutomaticData

isMC = False
year = '2016'
pol = 'MagDown'
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
        # Tau isolation BDT
		"BSTAUTAUTAUISOBDTFIRSTVALUETAUP" : "RELINFO('/Event/{0}/Phys/{1}/TauIsolationBDT', 'BSTAUTAUTAUISOBDTFIRSTVALUETAUP', -100)".format(stream,line),
		"BSTAUTAUTAUISOBDTSECONDVALUETAUP" : "RELINFO('/Event/{0}/Phys/{1}/TauIsolationBDT', 'BSTAUTAUTAUISOBDTSECONDVALUETAUP', -100)".format(stream,line),
		"BSTAUTAUTAUISOBDTTHIRDVALUETAUP" : "RELINFO('/Event/{0}/Phys/{1}/TauIsolationBDT', 'BSTAUTAUTAUISOBDTTHIRDVALUETAUP', -100)".format(stream,line),
		"BSTAUTAUTAUISOBDTFIRSTVALUETAUM" : "RELINFO('/Event/{0}/Phys/{1}/TauIsolationBDT', 'BSTAUTAUTAUISOBDTFIRSTVALUETAUM', -100)".format(stream,line),
		"BSTAUTAUTAUISOBDTSECONDVALUETAUM" : "RELINFO('/Event/{0}/Phys/{1}/TauIsolationBDT', 'BSTAUTAUTAUISOBDTSECONDVALUETAUM', -100)".format(stream,line),
		"BSTAUTAUTAUISOBDTTHIRDVALUETAUM" : "RELINFO('/Event/{0}/Phys/{1}/TauIsolationBDT', 'BSTAUTAUTAUISOBDTTHIRDVALUETAUM', -100)".format(stream,line),

		# Tau isolation
		"BSTAUTAUTAUISOFIRSTVALUETAUP" : "RELINFO('/Event/{0}/Phys/{1}/TauIsolation', 'BSTAUTAUTAUISOFIRSTVALUETAUP', -100)".format(stream,line),
		"BSTAUTAUTAUISOSECONDVALUETAUP" : "RELINFO('/Event/{0}/Phys/{1}/TauIsolation', 'BSTAUTAUTAUISOSECONDVALUETAUP', -100)".format(stream,line),
		"BSTAUTAUTAUISOFIRSTVALUETAUM" : "RELINFO('/Event/{0}/Phys/{1}/TauIsolation', 'BSTAUTAUTAUISOFIRSTVALUETAUM', -100)".format(stream,line),
		"BSTAUTAUTAUISOSECONDVALUETAUM" : "RELINFO('/Event/{0}/Phys/{1}/TauIsolation', 'BSTAUTAUTAUISOSECONDVALUETAUM', -100)".format(stream,line),

		# Track isolation BDT
		"BSTAUTAUTRKISOBDTFIRSTVALUETAUPPIM" : "RELINFO('/Event/{0}/Phys/{1}/TrackIsolationBDT', 'BSTAUTAUTRKISOBDTFIRSTVALUETAUPPIM', -100)".format(stream,line),
		"BSTAUTAUTRKISOBDTSECONDVALUETAUPPIM" : "RELINFO('/Event/{0}/Phys/{1}/TrackIsolationBDT', 'BSTAUTAUTRKISOBDTSECONDVALUETAUPPIM', -100)".format(stream,line),
		"BSTAUTAUTRKISOBDTTHIRDVALUETAUPPIM" : "RELINFO('/Event/{0}/Phys/{1}/TrackIsolationBDT', 'BSTAUTAUTRKISOBDTTHIRDVALUETAUPPIM', -100)".format(stream,line),
		"BSTAUTAUTRKISOBDTFIRSTVALUETAUPPIP1" : "RELINFO('/Event/{0}/Phys/{1}/TrackIsolationBDT', 'BSTAUTAUTRKISOBDTFIRSTVALUETAUPPIP1', -100)".format(stream,line),
		"BSTAUTAUTRKISOBDTSECONDVALUETAUPPIP1" : "RELINFO('/Event/{0}/Phys/{1}/TrackIsolationBDT', 'BSTAUTAUTRKISOBDTSECONDVALUETAUPPIP1', -100)".format(stream,line),
		"BSTAUTAUTRKISOBDTTHIRDVALUETAUPPIP1" : "RELINFO('/Event/{0}/Phys/{1}/TrackIsolationBDT', 'BSTAUTAUTRKISOBDTTHIRDVALUETAUPPIP1', -100)".format(stream,line),
		"BSTAUTAUTRKISOBDTFIRSTVALUETAUPPIP2" : "RELINFO('/Event/{0}/Phys/{1}/TrackIsolationBDT', 'BSTAUTAUTRKISOBDTFIRSTVALUETAUPPIP2', -100)".format(stream,line),
		"BSTAUTAUTRKISOBDTSECONDVALUETAUPPIP2" : "RELINFO('/Event/{0}/Phys/{1}/TrackIsolationBDT', 'BSTAUTAUTRKISOBDTSECONDVALUETAUPPIP2', -100)".format(stream,line),
		"BSTAUTAUTRKISOBDTTHIRDVALUETAUPPIP2" : "RELINFO('/Event/{0}/Phys/{1}/TrackIsolationBDT', 'BSTAUTAUTRKISOBDTTHIRDVALUETAUPPIP2', -100)".format(stream,line),
		"BSTAUTAUTRKISOBDTFIRSTVALUETAUMPIP" : "RELINFO('/Event/{0}/Phys/{1}/TrackIsolationBDT', 'BSTAUTAUTRKISOBDTFIRSTVALUETAUMPIP', -100)".format(stream,line),
		"BSTAUTAUTRKISOBDTSECONDVALUETAUMPIP" : "RELINFO('/Event/{0}/Phys/{1}/TrackIsolationBDT', 'BSTAUTAUTRKISOBDTSECONDVALUETAUMPIP', -100)".format(stream,line),
		"BSTAUTAUTRKISOBDTTHIRDVALUETAUMPIP" : "RELINFO('/Event/{0}/Phys/{1}/TrackIsolationBDT', 'BSTAUTAUTRKISOBDTTHIRDVALUETAUMPIP', -100)".format(stream,line),
		"BSTAUTAUTRKISOBDTFIRSTVALUETAUMPIM1" : "RELINFO('/Event/{0}/Phys/{1}/TrackIsolationBDT', 'BSTAUTAUTRKISOBDTFIRSTVALUETAUMPIM1', -100)".format(stream,line),
		"BSTAUTAUTRKISOBDTSECONDVALUETAUMPIM1" : "RELINFO('/Event/{0}/Phys/{1}/TrackIsolationBDT', 'BSTAUTAUTRKISOBDTSECONDVALUETAUMPIM1', -100)".format(stream,line),
		"BSTAUTAUTRKISOBDTTHIRDVALUETAUMPIM1" : "RELINFO('/Event/{0}/Phys/{1}/TrackIsolationBDT', 'BSTAUTAUTRKISOBDTTHIRDVALUETAUMPIM1', -100)".format(stream,line),
		"BSTAUTAUTRKISOBDTFIRSTVALUETAUMPIM2" : "RELINFO('/Event/{0}/Phys/{1}/TrackIsolationBDT', 'BSTAUTAUTRKISOBDTFIRSTVALUETAUMPIM2', -100)".format(stream,line),
		"BSTAUTAUTRKISOBDTSECONDVALUETAUMPIM2" : "RELINFO('/Event/{0}/Phys/{1}/TrackIsolationBDT', 'BSTAUTAUTRKISOBDTSECONDVALUETAUMPIM2', -100)".format(stream,line),
		"BSTAUTAUTRKISOBDTTHIRDVALUETAUMPIM2" : "RELINFO('/Event/{0}/Phys/{1}/TrackIsolationBDT', 'BSTAUTAUTRKISOBDTTHIRDVALUETAUMPIM2', -100)".format(stream,line),

		# Track isolation
		"BSTAUTAUTRKISOFIRSTVALUETAUPPIM" : "RELINFO('/Event/{0}/Phys/{1}/TrackIsolation', 'BSTAUTAUTRKISOFIRSTVALUETAUPPIM', -100)".format(stream,line),
		"BSTAUTAUTRKISOFIRSTVALUETAUPPIP1" : "RELINFO('/Event/{0}/Phys/{1}/TrackIsolation', 'BSTAUTAUTRKISOFIRSTVALUETAUPPIP1', -100)".format(stream,line),
		"BSTAUTAUTRKISOFIRSTVALUETAUPPIP2" : "RELINFO('/Event/{0}/Phys/{1}/TrackIsolation', 'BSTAUTAUTRKISOFIRSTVALUETAUPPIP2', -100)".format(stream,line),
		"BSTAUTAUTRKISOFIRSTVALUETAUMPIP" : "RELINFO('/Event/{0}/Phys/{1}/TrackIsolation', 'BSTAUTAUTRKISOFIRSTVALUETAUMPIP', -100)".format(stream,line),
		"BSTAUTAUTRKISOFIRSTVALUETAUMPIM1" : "RELINFO('/Event/{0}/Phys/{1}/TrackIsolation', 'BSTAUTAUTRKISOFIRSTVALUETAUMPIM1', -100)".format(stream,line),
		"BSTAUTAUTRKISOFIRSTVALUETAUMPIM2" : "RELINFO('/Event/{0}/Phys/{1}/TrackIsolation', 'BSTAUTAUTRKISOFIRSTVALUETAUMPIM2', -100)".format(stream,line),

		# CDFIso
		"BSTAUTAUCDFISO" : "RELINFO('/Event/{0}/Phys/{1}/CDFIso', 'BSTAUTAUCDFISO', -100)".format(stream,line),

		# ZVisoBDT
		"ZVISOTAUP" : "RELINFO('/Event/{0}/Phys/{1}/ZVisoBDT', 'ZVISOTAUP', -100)".format(stream,line),
		"ZVISOTAUM" : "RELINFO('/Event/{0}/Phys/{1}/ZVisoBDT', 'ZVISOTAUM', -100)".format(stream,line)
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
	# my decay fit (inputs for DECAY_FIT.C)
	# df = branch.addTupleTool("TupleToolDecayFit_DDK/df")
	# df.constrainToOriginVertex = True

    # standard DTF
    dtf = branch.addTupleTool("TupleToolDecayTreeFitter/dtf")
    dtf.daughtersToConstrain = ['D0', 'D_s+']
    dtf.UseFullTreeInName = True
    dtf.UpdateDaughters = True
    dtf.constrainToOriginVertex = True
    dtf.Substitutions = { 
        'B+ -> D~0 ^D+': 'D_s+', 
        'B- -> D0  ^D-': 'D_s-' 
    }

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
    addTrkInfo(DecayTreeTuple.D0bar_pi)
    addTrkInfo(DecayTreeTuple.Dsp_K1)
    addTrkInfo(DecayTreeTuple.Dsp_K2)
    addTrkInfo(DecayTreeTuple.Dsp_pi)
    #add M12, M23, M13 for tau's
    addMInfo(DecayTreeTuple.D0bar)
    addMInfo(DecayTreeTuple.Dsp)
    # add trigger information
    addTISTOS(DecayTreeTuple.Bp)
    if isMC:
        addTupleToolL0Calo(DecayTreeTuple.D0bar_K)
        addTupleToolL0Calo(DecayTreeTuple.D0bar_pi)
        addTupleToolL0Calo(DecayTreeTuple.Dsp_K1)
        addTupleToolL0Calo(DecayTreeTuple.Dsp_K2)
        addTupleToolL0Calo(DecayTreeTuple.Dsp_pi)

# Configure DaVinci
if isMC:
	DaVinci().TupleFile = 'DVntuple_MC_'+year+'_'+pol+'.root'
	# DaVinci().TupleFile = '/panfs/felician/BuD0Dps/DVntuple_MC_'+year+'_'+pol+'.root'
else:
	DaVinci().TupleFile = 'DVntuple_data_'+year+'_'+pol+'.root'
	# DaVinci().TupleFile = '/panfs/felician/BuD0Dps/DVntuple_data_'+year+'_'+pol+'.root'

# Stream and stripping line we want to use
if isMC:
	stream = 'AllStreams' # for filtered MC
else:	
	stream = 'Bhadron'
line = 'B2XTau_DD0_Line'

if isMC:
	DaVinci().Simulation = True

	# CondDBtag / DDDBtag specify the exact detector conditions with which the MC was generated. They are specified in the downloaded DST file.
	if year == '2016':
		DaVinci().DDDBtag = "dddb-20170721-3" 
		if pol == 'MagUp':
			DaVinci().CondDBtag = "sim-20170721-2-vc-mu100"  
		elif pol == 'MagDown':
			DaVinci().CondDBtag = "sim-20170721-2-vc-md100"
	elif year == '2017':
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
		b_strip    = AutomaticData('/Phys/{0}/Particles'.format(line)) 

	b_Filter    = FilterSelection("b_Filter", b_strip, Code = filtercode)
	dtt.Inputs = [b_Filter.outputLocation()]
else:
		dtt.Inputs = ['/Phys/{0}/Particles'.format(line)]

dtt.setDescriptorTemplate('${Bp}[B+ -> ${D0bar}(D~0 -> ${D0bar_K}K+ ${D0bar_pi}pi-) ${Dsp}(D+  -> ${Dsp_K1}K+ ${Dsp_K2}K- ${Dsp_pi}pi+)]CC')

# Generated (Acceptance cut: DaughtersInLHCb)
mcdtt.Decay = '[B+ => ^(D~0 => ^K+ ^pi-)  ^(D_s+ => ^K+ ^K- ^pi+)]CC' 

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
    STRIP_Code = "HLT_PASS_RE('StrippingB2XTau_DD0_LineDecision')"
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
# from GaudiConf import IOHelper
# IOHelper('ROOT').inputFiles([
# 	'/panfs/felician/BuD0D0K/DST_data_2018_MagUp/00076476_00010124_1.bhadron.mdst'
# 	# '/panfs/felician/BuD0Dps/DST_MC_2016/00125120_00000036_7.AllStreams.mdst',
# 	# '/panfs/felician/BuD0Dps/DST_MC_2016/00125120_00000107_7.AllStreams.mdst'
# ], clear=True)