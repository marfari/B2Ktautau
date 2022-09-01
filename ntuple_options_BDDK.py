from Configurables import DecayTreeTuple
from DecayTreeTuple.Configuration import *

isMC = False
year = '2016'
pol = 'MagUp'

# Stream and stripping line we want to use
if isMC:
	stream = 'AllStreams'
else:	
	stream = 'Bhadron'
line = 'B2DDKBeauty2CharmLine'

# Create an ntuple to capture B+ decays from the StrippingLine line
dtt = DecayTreeTuple('ntuple')

if isMC:
	dtt.Inputs = ['/Event/{0}/Phys/{1}/Particles'.format(stream, line)] # DST
else:
	dtt.Inputs = ['/Phys/{0}/Particles'.format(line)] # microDST 

dtt.Decay = '[B+ -> ^K+ ^(D+ -> ^K- ^pi+ ^pi+) ^(D- -> ^K+ ^pi- ^pi-) ]CC' 
dtt.addBranches({'Bp'  : '[B+ -> K+ (D+ -> K- pi+ pi+) (D- -> K+ pi- pi-) ]CC',
		'Kp'   : '[B+ -> ^K+ (D+ -> K- pi+ pi+) (D- -> K+ pi- pi-) ]CC', 
		'Dp' : '[B+ -> K+ ^(D+ -> K- pi+ pi+) (D- -> K+ pi- pi-) ]CC',
		'Dp_Km' : '[B+ -> K+ (D+ -> ^K- pi+ pi+) (D- -> K+ pi- pi-) ]CC',
		'Dp_pip1' : '[B+ -> K+ (D+ -> K- ^pi+ pi+) (D- -> K+ pi- pi-) ]CC',
		'Dp_pip2' : '[B+ -> K+ (D+ -> K- pi+ ^pi+) (D- -> K+ pi- pi-) ]CC',
		'Dm' : '[B+ -> K+ (D+ -> K- pi+ pi+) ^(D- -> K+ pi- pi-) ]CC',
		'Dm_Kp' : '[B+ -> K+ (D+ -> K- pi+ pi+) (D- -> ^K+ pi- pi-) ]CC',
		'Dm_pim1' : '[B+ -> K+ (D+ -> K- pi+ pi+) (D- -> K+ ^pi- pi-) ]CC',
		'Dm_pim2' : '[B+ -> K+ (D+ -> K- pi+ pi+) (D- -> K+ pi- ^pi-) ]CC'})

# add track info to particles which leave a track
dtt.addTupleTool('TupleToolTrackInfo')

if isMC:
	dtt.addTupleTool('TupleToolMCTruth')

# DecayTreeFitter
dtt.Bp.addTupleTool('TupleToolDecayTreeFitter/ConsBp')
dtt.Bp.ConsBp.Verbose = True
dtt.Dstar.ConsBp.daughtersToConstrain = ['D+', 'D-']
dtt.Bp.ConsBp.constrainToOriginVertex = True
dtt.Bp.ConsBp.UpdateDaughters = True

# Isolation information (stripping line output isolation variables)
from Configurables import LoKi__Hybrid__TupleTool
LoKi_Cone =  LoKi__Hybrid__TupleTool('LoKi_Cone')
LoKi_Cone.Variables = {
    "Bp_CONEANGLE_var1" : "RELINFO('/Phys/B2DDKBeauty2CharmLine/P2ConeVar1', 'CONEANGLE', -100)",
	"Bp_CONEAMULT_var1" : "RELINFO('/Phys/B2DDKBeauty2CharmLine/P2ConeVar1', 'CONEMULT', -100)",
	"Bp_CONEPTASYM_var1" : "RELINFO('/Phys/B2DDKBeauty2CharmLine/P2ConeVar1', 'CONEPTASYM', -100)",
	"Bp_CONEANGLE_var2" : "RELINFO('/Phys/B2DDKBeauty2CharmLine/P2ConeVar2', 'CONEANGLE', -100)",
	"Bp_CONEAMULT_var2" : "RELINFO('/Phys/B2DDKBeauty2CharmLine/P2ConeVar2', 'CONEMULT', -100)",
	"Bp_CONEPTASYM_var2" : "RELINFO('/Phys/B2DDKBeauty2CharmLine/P2ConeVar2', 'CONEPTASYM', -100)",
	"Bp_CONEANGLE_var3" : "RELINFO('/Phys/B2DDKBeauty2CharmLine/P2ConeVar3', 'CONEANGLE', -100)",
	"Bp_CONEAMULT_var3" : "RELINFO('/Phys/B2DDKBeauty2CharmLine/P2ConeVar3', 'CONEMULT', -100)",
	"Bp_CONEPTASYM_var3" : "RELINFO('/Phys/B2DDKBeauty2CharmLine/P2ConeVar3', 'CONEPTASYM', -100)"
}
dtt.Bp.addTupleTool('LoKi_Cone')

from Configurables import DaVinci

# Stripping filter
from PhysConf.Filters import LoKi_Filters
fltrs = LoKi_Filters (
    STRIP_Code = "HLT_PASS_RE('StrippingB2DDKBeauty2CharmLineDecision')"
)
DaVinci().EventPreFilters = fltrs.filters('Filters')

# Configure DaVinci
DaVinci().UserAlgorithms += [dtt]
if isMC:
	DaVinci().TupleFile = '/panfs/felician/B2DDK/ROOT_Sim/'+year+'/DVntuple_MC_'+year+'_'+pol+'.root'
	#DaVinci().TupleFile = 'DVntuple_MC_'+year+'_'+pol+'.root'
else:
	DaVinci().TupleFile = 'DVntuple_data_'+year+'_'+pol+'.root'
DaVinci().PrintFreq = 1000
DaVinci().DataType = year 
if isMC:
	DaVinci().Simulation = True
	DaVinci().InputType = 'DST'
	# CondDBtag / DDDBtag specify the exact detector conditions with which the MC was generated. They are specified in the downloaded DST file.
	if pol == 'MagUp':
		DaVinci().CondDBtag = "sim-20170721-2-vc-mu100"
	else:
		DaVinci().CondDBtag = "sim-20170721-2-vc-md100"
	DaVinci().DDDBtag = "dddb-20170721-3"
else:
	DaVinci().Simulation = False
	DaVinci().InputType = 'MDST'
	DaVinci().RootInTES = '/Event/{0}'.format(stream) # for MDST only
DaVinci().Lumi = not DaVinci().Simulation # Only ask for luminosity information when not using simulated data
	
DaVinci().EvtMax = -1 # how many events to run over (-1 = runs over all events)

from GaudiConf import IOHelper
# Use the local input data
#IOHelper().inputFiles([''''/panfs/felician/B2Ktautau/DST_Data/MagUp/2018/00092561_00000274_1.bhadron.mdst'], clear=True)

