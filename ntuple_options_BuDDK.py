from Configurables import DecayTreeTuple
from DecayTreeTuple.Configuration import *
from Configurables import LoKi__Hybrid__TupleTool
from Configurables import CheckPV
from Configurables import DaVinci

isMC = False
year = '2016'
pol = 'MagUp'

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
def addIsoInfo(branch, stream, line):
	LoKi_Cone =  branch.Bp.addTupleTool("LoKi::Hybrid::TupleTool/KoKi_Cone")

	LoKi_Cone.Variables = {
		"CONEANGLE_var1" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar1', 'CONEANGLE', -100)".format(stream,line),
		"CONEAMULT_var1" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar1', 'CONEMULT', -100)".format(stream,line),
		"CONEPTASYM_var1" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar1', 'CONEPTASYM', -100)".format(stream,line),
		"CONEANGLE_var2" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar2', 'CONEANGLE', -100)".format(stream,line),
		"CONEAMULT_var2" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar2', 'CONEMULT', -100)".format(stream,line),
		"CONEPTASYM_var2" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar2', 'CONEPTASYM', -100)".format(stream,line),
		"CONEANGLE_var3" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar3', 'CONEANGLE', -100)".format(stream,line),
		"CONEAMULT_var3" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar3', 'CONEMULT', -100)".format(stream,line),
		"CONEPTASYM_var3" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar3', 'CONEPTASYM', -100)".format(stream,line)
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
        "AM" : "PFUNA(AM)"
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
	dtt.Bp.addTupleTool('TupleToolDecayTreeFitter/ConsBp')
	dtt.Bp.ConsBp.Verbose = True
	dtt.Bp.ConsBp.daughtersToConstrain = ['D+', 'D-']
	dtt.Bp.ConsBp.constrainToOriginVertex = True
	dtt.Bp.ConsBp.UpdateDaughters = True

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
		"L0CALODecision"
		# HLT1
		"Hlt1TrackMVADecision",
		"Hlt1TrackMVALooseDecision",
		"Hlt1TwoTrackMVADecision",
		"Hlt1TwoTrackMVALooseDecision",
		"Hlt1CalibHighPTLowMultTrksDecision",
		"Hlt1IncPhiDecision",
		# HLT2
		"Hlt2ForwardDecision",
		"Hlt2XcMuXForTauB2XcFakeMuDecision",
		"Hlt2Topo2BodyDecision",
		"Hlt2Topo3BodyDecision",
		"Hlt2Topo4BodyDecision"
    ]

def addTools( DecayTreeTuple, stream, line ):
    #run modified DTF
    addKTTDTF(DecayTreeTuple.Bp) 
    #add some helpful variables for B0
    addTopInfo(DecayTreeTuple.Bp)
    #add isolation info
    addIsoInfo(DecayTreeTuple, stream, line)
    #add submasses
    addSubmassInfo(DecayTreeTuple.Bp)
    #add helpful info for final state tracks
    addTrkInfo(DecayTreeTuple.Dp_K)
    addTrkInfo(DecayTreeTuple.Dp_pi1)
    addTrkInfo(DecayTreeTuple.Dp_pi2)
    addTrkInfo(DecayTreeTuple.Dm_K)
    addTrkInfo(DecayTreeTuple.Dm_pi1)
    addTrkInfo(DecayTreeTuple.Dm_pi2)
    addTrkInfo(DecayTreeTuple.Kp)
    #add M12, M23, M13 for tau's
    addMInfo(DecayTreeTuple.Dp)
    addMInfo(DecayTreeTuple.Dm)
    # add trigger information
    addTISTOS(DecayTreeTuple.Bp)

# Stream and stripping line we want to use
if isMC:
	stream = 'AllStreams'
else:	
	stream = 'Bhadron'
line = 'B2DDKBeauty2CharmLine'
line_WS = 'B2DDKWSBeauty2CharmLine'

# Create an ntuple to capture B+ decays from the StrippingLine line
dtt = DecayTreeTuple('ntuple')
dtt_WS = DecayTreeTuple('ntuple_WS')
mcdtt = MCDecayTreeTuple('mc_ntuple')

if isMC:
	dtt.Inputs = ['/Event/{0}/Phys/{1}/Particles'.format(stream, line)] # DST
else:
	dtt.Inputs = ['/Phys/{0}/Particles'.format(line)] # microDST 
	dtt_WS.Inputs = ['/Phys/{0}/Particles'.format(line_WS)] # microDST 

dtt.setDescriptorTemplate('${Bp}[B+ -> ${Kp}K+ ${Dp}(D+ -> ${Dp_K}K- ${Dp_pi1}pi+ ${Dp_pi2}pi+) ${Dm}(D- -> ${Dm_K}K+ ${Dm_pi1}pi- ${Dm_pi2}pi-)]CC')
dtt_WS.setDescriptorTemplate('${Bp}[B+ -> ${Kp}K+ ${Dp}(D+ -> ${Dp_K}K- ${Dp_pi1}pi+ ${Dp_pi2}pi+) ${Dm}(D+ -> ${Dm_K}K- ${Dm_pi1}pi- ${Dm_pi2}pi-)]CC')
# Generated (Acceptance cut: DaughtersInLHCb)
mcdtt.Decay = '[B+ => ^(D+ => ^K- ^pi+ ^pi+)  ^(D- => ^K+ ^pi- ^pi-) K+]CC' 

# TUPLE TOOLS
ToolList = ['TupleToolTrackInfo', 'TupleToolVtxIsoln', 'TupleToolPrimaries', 'TupleToolPropertime', 'TupleToolAngles', 'TupleToolL0Calo', 'TupleToolRecoStats']
dtt.ToolList += ToolList
if not isMC:
	dtt_WS.ToolList += ToolList
if isMC:
	dtt.addTupleTool('TupleToolMCTruth')
	dtt.addTupleTool('TupleToolMCBackgroundInfo')

# Adds several tupletools, including DTF
addTools(dtt, stream, line)
if not isMC:
	addTools(dtt_WS, stream, line_WS)

GENToolList = ['MCTupleToolKinematic', 'MCTupleToolHierarchy', 'TupleToolEventInfo']
if isMC:
	mcdtt.ToolList += GENToolList

# Stripping filter
from PhysConf.Filters import LoKi_Filters
fltrs = LoKi_Filters (
    STRIP_Code = "HLT_PASS_RE('StrippingB2DDKBeauty2CharmLineDecision')"
)
if not isMC:
	DaVinci().EventPreFilters = fltrs.filters('Filters')

# PV checker -> TupleToolPrimaries requires there to be a PV in the event, otherwise it throws an error
DaVinci().UserAlgorithms += [CheckPV(), dtt]
if isMC:
	DaVinci().UserAlgorithms += [mcdtt]
if not isMC:
	DaVinci().UserAlgorithms += [dtt_WS]

# Configure DaVinci
if isMC:
	#DaVinci().TupleFile = '/panfs/felician/B2DDK/ROOT_Sim/'+year+'/DVntuple_MC_'+year+'_'+pol+'.root'
	DaVinci().TupleFile = 'DVntuple_MC_'+year+'_'+pol+'.root'
else:
	DaVinci().TupleFile = 'DVntuple_data_'+year+'_'+pol+'.root'
DaVinci().PrintFreq = 1000
DaVinci().DataType = year 
DaVinci().InputType = 'MDST'
if isMC:
	DaVinci().Simulation = True
	#DaVinci().InputType = 'DST'
	# CondDBtag / DDDBtag specify the exact detector conditions with which the MC was generated. They are specified in the downloaded DST file.
	if pol == 'MagUp':
		DaVinci().CondDBtag = "sim-20170721-2-vc-mu100"
	else:
		DaVinci().CondDBtag = "sim-20170721-2-vc-md100"
	DaVinci().DDDBtag = "dddb-20170721-3"
else:
	DaVinci().Simulation = False
	DaVinci().RootInTES = '/Event/{0}'.format(stream) # for MDST only
DaVinci().Lumi = not DaVinci().Simulation # Only ask for luminosity information when not using simulated data
	
DaVinci().EvtMax = 10000 # how many events to run over (-1 = runs over all events)

from GaudiConf import IOHelper
# Use the local input data
# IOHelper('ROOT').inputFiles([
# 'LFN:/lhcb/LHCb/Collision16/BHADRON.MDST/00070442/0000/00070442_00000166_1.bhadron.mdst']