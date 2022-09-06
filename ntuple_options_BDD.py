# ntuple options for B+ -> D0barD+ (D0bar -> a1+K- -> 3pi+ K-) (D+ -> a1+ K0 -> 3pi+ K0)
# decay is reconstructed as if it were our sginal B -> K tautau

from Configurables import DecayTreeTuple
from DecayTreeTuple.Configuration import *

isMC = True
year = '2016'
pol = 'MagDown'

# Stream and stripping line we want to use
if isMC:
	stream = 'AllStreams'
else:	
	stream = 'Bhadron'
line = 'B2KTauTauLine'

# Create an ntuple to capture B+ decays from the StrippingLine line
dtt = DecayTreeTuple('ntuple')
mcdtt = MCDecayTreeTuple('mc_ntuple')

if isMC:
	dtt.Inputs = ['/Event/{0}/Phys/{1}/Particles'.format(stream, line)] # DST
else:
	dtt.Inputs = ['/Phys/{0}/Particles'.format(line)] # microDST 

# Reconstructed
dtt.Decay = '[B+ -> ^K+ ^(tau+ -> ^pi+ ^pi- ^pi+) ^(tau- -> ^pi- ^pi+ ^pi-) ]CC' 
dtt.addBranches({'Bp'  : '[B+ -> K+ (tau+ -> pi+ pi- pi+) (tau- -> pi- pi+ pi-) ]CC',
		'Kp'   : '[B+ -> ^K+ (tau+ -> pi+ pi- pi+) (tau- -> pi- pi+ pi-) ]CC', 
		'taup' : '[B+ -> K+ ^(tau+ -> pi+ pi- pi+) (tau- -> pi- pi+ pi-) ]CC',
		'taup_pip0' : '[B+ -> K+ (tau+ -> ^pi+ pi- pi+) (tau- -> pi- pi+ pi-) ]CC',
		'taup_pim0' : '[B+ -> K+ (tau+ -> pi+ ^pi- pi+) (tau- -> pi- pi+ pi-) ]CC',
		'taup_pip1' : '[B+ -> K+ (tau+ -> pi+ pi- ^pi+) (tau- -> pi- pi+ pi-) ]CC',
		'taum' : '[B+ -> K+ (tau+ -> pi+ pi- pi+) ^(tau- -> pi- pi+ pi-) ]CC',
		'taum_pim0' : '[B+ -> K+ (tau+ -> pi+ pi- pi+) (tau- -> ^pi- pi+ pi-) ]CC',
		'taum_pip0' : '[B+ -> K+ (tau+ -> pi+ pi- pi+) (tau- -> pi- ^pi+ pi-) ]CC',
		'taum_pim1' : '[B+ -> K+ (tau+ -> pi+ pi- pi+) (tau- -> pi- pi+ ^pi-) ]CC'})
#dtt.setDescriptorTemplate('${B+}[B+ -> ${K+}K+ ${tau+}(tau+ -> ${pi+}pi+ ${pi-}pi- ${pi+}pi+) ${tau-}(tau- -> ${pi-}pi- ${pi+}pi+ ${pi-}pi-) ]CC')

# Generated (Acceptance cut: DaughtersInLHCb)
mcdtt.Decay = '[B+ => ^(D~0 => ^(a_1(1260)- => ^(rho(770)0 => ^pi+ ^pi-) ^pi-) ^K+)  ^(D+ => ^(a_1(1260)+ => ^(rho(770)0 => ^pi+ ^pi-) ^pi+) ^K0) ]CC' 

# add track info to particles which leave a track
dtt.addTupleTool('TupleToolTrackInfo')

if isMC:
	dtt.addTupleTool('TupleToolMCTruth')

mcdtt.ToolList = []
mcdtt.addTupleTool('LoKi::Hybrid::MCTupleTool/basicLoKiTT')

# DecayTreeFitter
dtt.Bp.addTupleTool('B2KtautauDTF1/ConsBp')
dtt.Bp.ConsBp.Verbose = True
dtt.Bp.ConsBp.constrainToOriginVertex = True
dtt.Bp.ConsBp.daughtersToConstrain = ['tau+', 'tau-', "nu_tau", "nu_tau~"]
dtt.Bp.ConsBp.UpdateDaughters = True
#B+ is constrained to come from PV + tau+ and tau- are constrained to have the tau PDG mass

from Configurables import DaVinci

# from PhysConf.Filters import LoKi_Filters
# fltrs = LoKi_Filters (
#     STRIP_Code = "HLT_PASS_RE('StrippingB2KTauTauLineDecision')"
# )
# DaVinci().EventPreFilters = fltrs.filters('Filters')

# Configure DaVinci
DaVinci().UserAlgorithms += [dtt]
DaVinci().UserAlgorithms += [mcdtt]
if isMC:
	DaVinci().TupleFile = '/panfs/felician/BDD_bkg/ROOT_Sim/'+year+'/DVntuple_MC_'+year+'_'+pol+'.root'
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

IOHelper().inputFiles(['/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_902135752.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_902135759.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_902135762.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_902135763.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_902135766.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_902135767.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_902135768.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_902135769.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_902135770.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_902135771.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_902135772.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_902135774.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_902135775.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_902135776.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_902135777.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_902135779.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_902135781.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_902135782.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_902135783.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_902135784.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_902135786.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_902135788.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_902135789.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_902135790.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_902135791.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_902135792.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_902135793.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_902135796.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_902135798.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_902135799.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170417.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170419.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170427.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170436.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170453.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170519.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170521.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170529.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170542.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170558.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170560.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170588.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170612.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170622.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170639.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170679.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170725.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170747.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170760.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170783.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170836.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170839.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170842.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170860.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170866.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170882.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170898.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170910.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170912.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170917.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170944.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830170964.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830171004.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830171009.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830171025.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830171029.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_830171035.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_902135792.dst'])