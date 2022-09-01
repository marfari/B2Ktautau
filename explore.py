import sys

import GaudiPython as GP
from GaudiConf import IOHelper
from Configurables import DaVinci

dv = DaVinci()
dv.DataType = '2016'
dv.Simulation = True

# Retrieve file path to open as the last command line argument
#inputFiles = [sys.argv[-1]] # takes it from command line

inputFiles = [ "/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172566.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172567.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172569.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172569.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172570.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172571.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172572.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172573.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172574.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172575.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172576.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172577.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172578.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172579.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172580.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172581.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172582.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172583.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172584.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172585.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172586.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172587.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172588.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172589.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172590.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172591.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172592.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172593.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172594.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172595.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172596.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172597.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172598.dst",
"/panfs/avenkate/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_413172599.dst"
]

IOHelper('ROOT').inputFiles(inputFiles)

appMgr = GP.AppMgr()
evt = appMgr.evtsvc()

appMgr.run(1)

def nodes(evt, root='/Event'):
    """List all nodes in `evt` starting from `node`."""
    node_list = [root]
    for leaf in evt.leaves(evt[root]):
        node_list += nodes(evt, leaf.identifier())
    return node_list

def advance(decision):
    """Advance until stripping decision is true, returns
    number of events by which we advanced"""
    n = 0
    while evt['/Event']:
        reports = evt['/Event/Strip/Phys/DecReports']
        if reports.hasDecisionName('Stripping{}Decision'.format(decision)):
            break

        appMgr.run(1)
        n += 1
    return n







