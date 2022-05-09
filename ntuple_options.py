from Configurables import DecayTreeTuple
from DecayTreeTuple.Configuration import *

# Stream and stripping line we want to use
stream = 'Bhadron'
line = 'B2KTauTauLine'

# Create an ntuple to capture B+ decays from the StrippingLine line
dtt = DecayTreeTuple('ntuple')
dtt.Inputs = ['/Phys/{0}/Particles'.format(line)] # change for microDST / DST

dtt.Decay = '[B+ -> ^K+ ^(tau+ -> ^pi+ ^pi- ^pi+) ^(tau- -> ^pi- ^pi+ ^pi-) ]CC' 
dtt.addBranches({'Bp'  : '[B+ -> K+ (tau+ -> pi+ pi- pi+) (tau- -> pi- pi+ pi-) ]CC',
		'Kp'   : '[B+ -> ^K+ (tau+ -> pi+ pi- pi+) (tau- -> pi- pi+ pi-) ]CC', 
		'taup' : '[B+ -> K+ ^(tau+ -> pi+ pi- pi+) (tau- -> pi- pi+ pi-) ]CC',
		'pip0' : '[B+ -> K+ (tau+ -> ^pi+ pi- pi+) (tau- -> pi- pi+ pi-) ]CC',
		'pim0' : '[B+ -> K+ (tau+ -> pi+ ^pi- pi+) (tau- -> pi- pi+ pi-) ]CC',
		'pip1' : '[B+ -> K+ (tau+ -> pi+ pi- ^pi+) (tau- -> pi- pi+ pi-) ]CC',
		'taum' : '[B+ -> K+ (tau+ -> pi+ pi- pi+) ^(tau- -> pi- pi+ pi-) ]CC',
		'pim1' : '[B+ -> K+ (tau+ -> pi+ pi- pi+) (tau- -> ^pi- pi+ pi-) ]CC',
		'pip2' : '[B+ -> K+ (tau+ -> pi+ pi- pi+) (tau- -> pi- ^pi+ pi-) ]CC',
		'pim2' : '[B+ -> K+ (tau+ -> pi+ pi- pi+) (tau- -> pi- pi+ ^pi-) ]CC'})
#dtt.setDescriptorTemplate('${B+}[B+ -> ${K+}K+ ${tau+}(tau+ -> ${pi+}pi+ ${pi-}pi- ${pi+}pi+) ${tau-}(tau- -> ${pi-}pi- ${pi+}pi+ ${pi-}pi-) ]CC')

# add track info to particles which leave a track
dtt.addTupleTool('TupleToolTrackInfo')

# DecayTreeFitter
dtt.Bp.addTupleTool('TupleToolDecayTreeFitter/ConsBp')
dtt.Bp.ConsBp.Verbose = True
dtt.Bp.ConsBp.constrainToOriginVertex = True
dtt.Bp.ConsBp.daughtersToConstrain = ['tau+', 'tau-']
dtt.Bp.ConsBp.UpdateDaughters = True
# B+ is constrained to come from PV + tau+ and tau- are constrained to have the tau PDG mass

# missing: neutrinos

from Configurables import DaVinci
DaVinci().RootInTES = '/Event/{0}'.format(stream) # for MDST only

from PhysConf.Filters import LoKi_Filters
fltrs = LoKi_Filters (
    STRIP_Code = "HLT_PASS_RE('StrippingB2KTauTauLineDecision')"
)
DaVinci().EventPreFilters = fltrs.filters('Filters')

# Configure DaVinci
DaVinci().UserAlgorithms += [dtt]
DaVinci().InputType = 'MDST'
DaVinci().TupleFile = 'DVntuple_data_2018.root'
DaVinci().PrintFreq = 100
DaVinci().DataType = '2018' # collision year (change)
DaVinci().Simulation = False # simulation or data? (change)
# Only ask for luminosity information when not using simulated data
DaVinci().Lumi = not DaVinci().Simulation
DaVinci().EvtMax = 1000 # how many events to run over (-1 = runs over all events)
# CondDBtag / DDDBtag specify the exatct detector conditions with which the MC was generates (change). They are specified in the downloaded DST file.
# Comment them for real data
#DaVinci().CondDBtag = 'cond-20180202'
#DaVinci().DDDBtag = 'dddb-20190206-3' 

from GaudiConf import IOHelper

# Use the local input data
#IOHelper().inputFiles(['/panfs/felician/B2Ktautau/DST_Data/MagDown/2018/00092561_00000274_1.bhadron.mdst'], clear=True)
