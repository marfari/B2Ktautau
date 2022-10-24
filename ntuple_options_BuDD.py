# ntuple options for B+ -> D0barD+ (D0bar -> a1+K- -> 3pi+ K-) (D+ -> a1+ K0 -> 3pi+ K0)
# decay is reconstructed as if it were our sginal B -> K tautau

from Configurables import DecayTreeTuple
from DecayTreeTuple.Configuration import *
from Configurables import LoKi__Hybrid__TupleTool
from Configurables import CheckPV
from Configurables import DaVinci
from PhysConf.Selections import FilterSelection
from PhysSelPython.Wrappers import AutomaticData

isMC = True
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
def addIsoInfo(branch, stream, line, year):
	LoKi_Cone =  branch.Bp.addTupleTool("LoKi::Hybrid::TupleTool/KoKi_Cone")

	if(year == '2018'):
		LoKi_Cone.Variables = {
			"VTXISONUMVTX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_VertexIsoInfo', 'VTXISONUMVTX', -100)".format(stream,line),
			"VTXISODCHI2ONETRACK_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_VertexIsoInfo', 'VTXISODCHI2ONETRACK', -100)".format(stream,line),
			"VTXISODCHI2MASSONETRACK_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_VertexIsoInfo', 'VTXISODCHI2MASSONETRACK', -100)".format(stream,line),
			"VTXISODCHI2TWOTRACK_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_VertexIsoInfo', 'VTXISODCHI2TWOTRACK', -100)".format(stream,line),
			"VTXISODCHI2MASSTWOTRACK_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_VertexIsoInfo', 'VTXISODCHI2MASSTWOTRACK', -100)".format(stream,line),

			"VTXISONUMVTX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_VertexIsoInfo', 'VTXISONUMVTX', -100)".format(stream,line),
			"VTXISODCHI2ONETRACK_K" : "RELINFO('/Event/{0}/Phys/{1}/H_VertexIsoInfo', 'VTXISODCHI2ONETRACK', -100)".format(stream,line),
			"VTXISODCHI2MASSONETRACK_K" : "RELINFO('/Event/{0}/Phys/{1}/H_VertexIsoInfo', 'VTXISODCHI2MASSONETRACK', -100)".format(stream,line),
			"VTXISODCHI2TWOTRACK_K" : "RELINFO('/Event/{0}/Phys/{1}/H_VertexIsoInfo', 'VTXISODCHI2TWOTRACK', -100)".format(stream,line),
			"VTXISODCHI2MASSTWOTRACK_K" : "RELINFO('/Event/{0}/Phys/{1}/H_VertexIsoInfo', 'VTXISODCHI2MASSTWOTRACK', -100)".format(stream,line),

			"VTXISONUMVTX_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_VertexIsoInfo', 'VTXISONUMVTX', -100)".format(stream,line),
			"VTXISODCHI2ONETRACK_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_VertexIsoInfo', 'VTXISODCHI2ONETRACK', -100)".format(stream,line),
			"VTXISODCHI2MASSONETRACK_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_VertexIsoInfo', 'VTXISODCHI2MASSONETRACK', -100)".format(stream,line),
			"VTXISODCHI2TWOTRACK_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_VertexIsoInfo', 'VTXISODCHI2TWOTRACK', -100)".format(stream,line),
			"VTXISODCHI2MASSTWOTRACK_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_VertexIsoInfo', 'VTXISODCHI2MASSTWOTRACK', -100)".format(stream,line),

			"VTXISONUMVTX_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_VertexIsoInfo', 'VTXISONUMVTX', -100)".format(stream,line),
			"VTXISODCHI2ONETRACK_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_VertexIsoInfo', 'VTXISODCHI2ONETRACK', -100)".format(stream,line),
			"VTXISODCHI2MASSONETRACK_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_VertexIsoInfo', 'VTXISODCHI2MASSONETRACK', -100)".format(stream,line),
			"VTXISODCHI2TWOTRACK_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_VertexIsoInfo', 'VTXISODCHI2TWOTRACK', -100)".format(stream,line),
			"VTXISODCHI2MASSTWOTRACK_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_VertexIsoInfo', 'VTXISODCHI2MASSTWOTRACK', -100)".format(stream,line),

			"CC_ANGLE_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'CC_ANGLE', -100.)".format(stream, line),
			"CC_MULT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'CC_MULT', -100.)".format(stream, line),
			"CC_SPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'CC_SPT', -100.)".format(stream, line),
			"CC_VPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'CC_VPT', -100.)".format(stream, line),
			"CC_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'CC_PX', -100.)".format(stream, line),
			"CC_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'CC_PY', -100.)".format(stream, line),
			"CC_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'CC_PZ', -100.)".format(stream, line),
			"CC_PASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'CC_PASYM', -100.)".format(stream, line),
			"CC_PTASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'CC_PTASYM', -100.)".format(stream, line),
			"CC_PXASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'CC_PXASYM', -100.)".format(stream, line),
			"CC_PYASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'CC_PYASYM', -100.)".format(stream, line),
			"CC_PZASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'CC_PZASYM', -100.)".format(stream, line),
			"CC_DELTAETA_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'CC_DELTAETA', -100.)".format(stream, line),
			"CC_DELTAPHI_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'CC_DELTAPHI', -100.)".format(stream, line),
			"CC_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'CC_PX', -100.)".format(stream, line),
			"CC_IT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'CC_IT', -100.)".format(stream, line),
			"CC_MAXPT_Q_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'CC_MAXPT_Q', -100.)".format(stream, line),
			"CC_MAXPT_PT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'CC_MAXPT_PT', -100.)".format(stream, line),
			"CC_MAXPT_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'CC_MAXPT_PX', -100.)".format(stream, line),
			"CC_MAXPT_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'CC_MAXPT_PY', -100.)".format(stream, line),
			"CC_MAXPT_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'CC_MAXPT_PZ', -100.)".format(stream, line),
			"CC_MAXPT_PE_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'CC_MAXPT_PE', -100.)".format(stream, line),
			"NC_ANGLE_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'NC_ANGLE', -100.)".format(stream, line),
			"NC_MULT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'NC_MULT', -100.)".format(stream, line),
			"NC_SPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'NC_SPT', -100.)".format(stream, line),
			"NC_VPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'NC_VPT', -100.)".format(stream, line),
			"NC_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'NC_PX', -100.)".format(stream, line),
			"NC_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'NC_PY', -100.)".format(stream, line),
			"NC_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'NC_PZ', -100.)".format(stream, line),
			"NC_PASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'NC_PASYM', -100.)".format(stream, line),
			"NC_PTASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'NC_PTASYM', -100.)".format(stream, line),
			"NC_PXASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'NC_PXASYM', -100.)".format(stream, line),
			"NC_PYASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'NC_PYASYM', -100.)".format(stream, line),
			"NC_PZASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'NC_PZASYM', -100.)".format(stream, line),
			"NC_DELTAETA_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'NC_DELTAETA', -100.)".format(stream, line),
			"NC_DELTAPHI_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'NC_DELTAPHI', -100.)".format(stream, line),
			"NC_IT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'NC_IT', -100.)".format(stream, line),
			"NC_MAXPT_PT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'NC_MAXPT_PT', -100.)".format(stream, line),
			"NC_MAXPT_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'NC_MAXPT_PX', -100.)".format(stream, line),
			"NC_MAXPT_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'NC_MAXPT_PY', -100.)".format(stream, line),
			"NC_MAXPT_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'NC_MAXPT_PZ', -100.)".format(stream, line),
			"CCNC_IT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo', 'CCNC_IT', -100.)".format(stream, line),

			"CC_ANGLE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'CC_ANGLE', -100.)".format(stream, line),
			"CC_MULT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'CC_MULT', -100.)".format(stream, line),
			"CC_SPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'CC_SPT', -100.)".format(stream, line),
			"CC_VPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'CC_VPT', -100.)".format(stream, line),
			"CC_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'CC_PX', -100.)".format(stream, line),
			"CC_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'CC_PY', -100.)".format(stream, line),
			"CC_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'CC_PZ', -100.)".format(stream, line),
			"CC_PASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'CC_PASYM', -100.)".format(stream, line),
			"CC_PTASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'CC_PTASYM', -100.)".format(stream, line),
			"CC_PXASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'CC_PXASYM', -100.)".format(stream, line),
			"CC_PYASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'CC_PYASYM', -100.)".format(stream, line),
			"CC_PZASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'CC_PZASYM', -100.)".format(stream, line),
			"CC_DELTAETA_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'CC_DELTAETA', -100.)".format(stream, line),
			"CC_DELTAPHI_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'CC_DELTAPHI', -100.)".format(stream, line),
			"CC_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'CC_PX', -100.)".format(stream, line),
			"CC_IT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'CC_IT', -100.)".format(stream, line),
			"CC_MAXPT_Q_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'CC_MAXPT_Q', -100.)".format(stream, line),
			"CC_MAXPT_PT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'CC_MAXPT_PT', -100.)".format(stream, line),
			"CC_MAXPT_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'CC_MAXPT_PX', -100.)".format(stream, line),
			"CC_MAXPT_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'CC_MAXPT_PY', -100.)".format(stream, line),
			"CC_MAXPT_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'CC_MAXPT_PZ', -100.)".format(stream, line),
			"CC_MAXPT_PE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'CC_MAXPT_PE', -100.)".format(stream, line),
			"NC_ANGLE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'NC_ANGLE', -100.)".format(stream, line),
			"NC_MULT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'NC_MULT', -100.)".format(stream, line),
			"NC_SPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'NC_SPT', -100.)".format(stream, line),
			"NC_VPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'NC_VPT', -100.)".format(stream, line),
			"NC_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'NC_PX', -100.)".format(stream, line),
			"NC_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'NC_PY', -100.)".format(stream, line),
			"NC_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'NC_PZ', -100.)".format(stream, line),
			"NC_PASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'NC_PASYM', -100.)".format(stream, line),
			"NC_PTASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'NC_PTASYM', -100.)".format(stream, line),
			"NC_PXASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'NC_PXASYM', -100.)".format(stream, line),
			"NC_PYASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'NC_PYASYM', -100.)".format(stream, line),
			"NC_PZASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'NC_PZASYM', -100.)".format(stream, line),
			"NC_DELTAETA_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'NC_DELTAETA', -100.)".format(stream, line),
			"NC_DELTAPHI_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'NC_DELTAPHI', -100.)".format(stream, line),
			"NC_IT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'NC_IT', -100.)".format(stream, line),
			"NC_MAXPT_PT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'NC_MAXPT_PT', -100.)".format(stream, line),
			"NC_MAXPT_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'NC_MAXPT_PX', -100.)".format(stream, line),
			"NC_MAXPT_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'NC_MAXPT_PY', -100.)".format(stream, line),
			"NC_MAXPT_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'NC_MAXPT_PZ', -100.)".format(stream, line),
			"CCNC_IT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo', 'CCNC_IT', -100.)".format(stream, line),
	
			"CC_ANGLE_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'CC_ANGLE', -100.)".format(stream, line),
			"CC_MULT_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'CC_MULT', -100.)".format(stream, line),
			"CC_SPT_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'CC_SPT', -100.)".format(stream, line),
			"CC_VPT_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'CC_VPT', -100.)".format(stream, line),
			"CC_PX_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'CC_PX', -100.)".format(stream, line),
			"CC_PY_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'CC_PY', -100.)".format(stream, line),
			"CC_PZ_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'CC_PZ', -100.)".format(stream, line),
			"CC_PASYM_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'CC_PASYM', -100.)".format(stream, line),
			"CC_PTASYM_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'CC_PTASYM', -100.)".format(stream, line),
			"CC_PXASYM_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'CC_PXASYM', -100.)".format(stream, line),
			"CC_PYASYM_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'CC_PYASYM', -100.)".format(stream, line),
			"CC_PZASYM_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'CC_PZASYM', -100.)".format(stream, line),
			"CC_DELTAETA_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'CC_DELTAETA', -100.)".format(stream, line),
			"CC_DELTAPHI_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'CC_DELTAPHI', -100.)".format(stream, line),
			"CC_PX_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'CC_PX', -100.)".format(stream, line),
			"CC_IT_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'CC_IT', -100.)".format(stream, line),
			"CC_MAXPT_Q_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'CC_MAXPT_Q', -100.)".format(stream, line),
			"CC_MAXPT_PT_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'CC_MAXPT_PT', -100.)".format(stream, line),
			"CC_MAXPT_PX_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'CC_MAXPT_PX', -100.)".format(stream, line),
			"CC_MAXPT_PY_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'CC_MAXPT_PY', -100.)".format(stream, line),
			"CC_MAXPT_PZ_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'CC_MAXPT_PZ', -100.)".format(stream, line),
			"CC_MAXPT_PE_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'CC_MAXPT_PE', -100.)".format(stream, line),
			"NC_ANGLE_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'NC_ANGLE', -100.)".format(stream, line),
			"NC_MULT_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'NC_MULT', -100.)".format(stream, line),
			"NC_SPT_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'NC_SPT', -100.)".format(stream, line),
			"NC_VPT_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'NC_VPT', -100.)".format(stream, line),
			"NC_PX_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'NC_PX', -100.)".format(stream, line),
			"NC_PY_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'NC_PY', -100.)".format(stream, line),
			"NC_PZ_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'NC_PZ', -100.)".format(stream, line),
			"NC_PASYM_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'NC_PASYM', -100.)".format(stream, line),
			"NC_PTASYM_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'NC_PTASYM', -100.)".format(stream, line),
			"NC_PXASYM_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'NC_PXASYM', -100.)".format(stream, line),
			"NC_PYASYM_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'NC_PYASYM', -100.)".format(stream, line),
			"NC_PZASYM_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'NC_PZASYM', -100.)".format(stream, line),
			"NC_DELTAETA_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'NC_DELTAETA', -100.)".format(stream, line),
			"NC_DELTAPHI_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'NC_DELTAPHI', -100.)".format(stream, line),
			"NC_IT_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'NC_IT', -100.)".format(stream, line),
			"NC_MAXPT_PT_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'NC_MAXPT_PT', -100.)".format(stream, line),
			"NC_MAXPT_PX_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'NC_MAXPT_PX', -100.)".format(stream, line),
			"NC_MAXPT_PY_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'NC_MAXPT_PY', -100.)".format(stream, line),
			"NC_MAXPT_PZ_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'NC_MAXPT_PZ', -100.)".format(stream, line),
			"CCNC_IT_tau1" : "RELINFO('/Event/{0}/Phys/{1}/Tau1_ConeIsoInfo', 'CCNC_IT', -100.)".format(stream, line),

			"CC_ANGLE_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'CC_ANGLE', -100.)".format(stream, line),
			"CC_MULT_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'CC_MULT', -100.)".format(stream, line),
			"CC_SPT_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'CC_SPT', -100.)".format(stream, line),
			"CC_VPT_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'CC_VPT', -100.)".format(stream, line),
			"CC_PX_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'CC_PX', -100.)".format(stream, line),
			"CC_PY_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'CC_PY', -100.)".format(stream, line),
			"CC_PZ_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'CC_PZ', -100.)".format(stream, line),
			"CC_PASYM_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'CC_PASYM', -100.)".format(stream, line),
			"CC_PTASYM_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'CC_PTASYM', -100.)".format(stream, line),
			"CC_PXASYM_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'CC_PXASYM', -100.)".format(stream, line),
			"CC_PYASYM_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'CC_PYASYM', -100.)".format(stream, line),
			"CC_PZASYM_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'CC_PZASYM', -100.)".format(stream, line),
			"CC_DELTAETA_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'CC_DELTAETA', -100.)".format(stream, line),
			"CC_DELTAPHI_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'CC_DELTAPHI', -100.)".format(stream, line),
			"CC_PX_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'CC_PX', -100.)".format(stream, line),
			"CC_IT_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'CC_IT', -100.)".format(stream, line),
			"CC_MAXPT_Q_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'CC_MAXPT_Q', -100.)".format(stream, line),
			"CC_MAXPT_PT_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'CC_MAXPT_PT', -100.)".format(stream, line),
			"CC_MAXPT_PX_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'CC_MAXPT_PX', -100.)".format(stream, line),
			"CC_MAXPT_PY_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'CC_MAXPT_PY', -100.)".format(stream, line),
			"CC_MAXPT_PZ_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'CC_MAXPT_PZ', -100.)".format(stream, line),
			"CC_MAXPT_PE_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'CC_MAXPT_PE', -100.)".format(stream, line),
			"NC_ANGLE_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'NC_ANGLE', -100.)".format(stream, line),
			"NC_MULT_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'NC_MULT', -100.)".format(stream, line),
			"NC_SPT_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'NC_SPT', -100.)".format(stream, line),
			"NC_VPT_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'NC_VPT', -100.)".format(stream, line),
			"NC_PX_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'NC_PX', -100.)".format(stream, line),
			"NC_PY_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'NC_PY', -100.)".format(stream, line),
			"NC_PZ_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'NC_PZ', -100.)".format(stream, line),
			"NC_PASYM_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'NC_PASYM', -100.)".format(stream, line),
			"NC_PTASYM_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'NC_PTASYM', -100.)".format(stream, line),
			"NC_PXASYM_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'NC_PXASYM', -100.)".format(stream, line),
			"NC_PYASYM_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'NC_PYASYM', -100.)".format(stream, line),
			"NC_PZASYM_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'NC_PZASYM', -100.)".format(stream, line),
			"NC_DELTAETA_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'NC_DELTAETA', -100.)".format(stream, line),
			"NC_DELTAPHI_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'NC_DELTAPHI', -100.)".format(stream, line),
			"NC_IT_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'NC_IT', -100.)".format(stream, line),
			"NC_MAXPT_PT_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'NC_MAXPT_PT', -100.)".format(stream, line),
			"NC_MAXPT_PX_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'NC_MAXPT_PX', -100.)".format(stream, line),
			"NC_MAXPT_PY_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'NC_MAXPT_PY', -100.)".format(stream, line),
			"NC_MAXPT_PZ_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'NC_MAXPT_PZ', -100.)".format(stream, line),
			"CCNC_IT_tau2" : "RELINFO('/Event/{0}/Phys/{1}/Tau2_ConeIsoInfo', 'CCNC_IT', -100.)".format(stream, line)
		}
	else:
		LoKi_Cone.Variables = {
			"VTXISONUMVTX_B" : "RELINFO('/Event/{0}/Phys/{1}/VertexIsoInfo', 'VTXISONUMVTX', -100)".format(stream,line),
			"VTXISODCHI2ONETRACK_B" : "RELINFO('/Event/{0}/Phys/{1}/VertexIsoInfo', 'VTXISODCHI2ONETRACK', -100)".format(stream,line),
			"VTXISODCHI2MASSONETRACK_B" : "RELINFO('/Event/{0}/Phys/{1}/VertexIsoInfo', 'VTXISODCHI2MASSONETRACK', -100)".format(stream,line),
			"VTXISODCHI2TWOTRACK_B" : "RELINFO('/Event/{0}/Phys/{1}/VertexIsoInfo', 'VTXISODCHI2TWOTRACK', -100)".format(stream,line),
			"VTXISODCHI2MASSTWOTRACK_B" : "RELINFO('/Event/{0}/Phys/{1}/VertexIsoInfo', 'VTXISODCHI2MASSTWOTRACK', -100)".format(stream,line),

			"CONEANGLE_B" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar05_B', 'CONEANGLE', -100.)".format(stream, line),
			"CONEMULT_B" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar05_B', 'CONEMULT', -100.)".format(stream, line),
			"CONEPASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar05_B', 'CONEPASYM', -100.)".format(stream, line),
			"CONEPTASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar05_B', 'CONEPTASYM', -100.)".format(stream, line),
			"CONEDELTAETA_B" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar05_B', 'CONEDELTAETA', -100.)".format(stream, line),

			"CONEANGLE_K" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar05_X', 'CONEANGLE', -100.)".format(stream, line),
			"CONEMULT_K" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar05_X', 'CONEMULT', -100.)".format(stream, line),
			"CONEPASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar05_X', 'CONEPASYM', -100.)".format(stream, line),
			"CONEPTASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar05_X', 'CONEPTASYM', -100.)".format(stream, line),
			"CONEDELTAETA_K" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar05_X', 'CONEDELTAETA', -100.)".format(stream, line),

			"CONEANGLE_tau1" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar05_Tau1', 'CONEANGLE', -100.)".format(stream, line),
			"CONEMULT_tau1" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar05_Tau1', 'CONEMULT', -100.)".format(stream, line),
			"CONEPASYM_tau1" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar05_Tau1', 'CONEPASYM', -100.)".format(stream, line),
			"CONEPTASYM_tau1" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar05_Tau1', 'CONEPTASYM', -100.)".format(stream, line),
			"CONEDELTAETA_tau1" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar05_Tau1', 'CONEDELTAETA', -100.)".format(stream, line),

			"CONEANGLE_tau2" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar05_Tau2', 'CONEANGLE', -100.)".format(stream, line),
			"CONEMULT_tau2" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar05_Tau2', 'CONEMULT', -100.)".format(stream, line),
			"CONEPASYM_tau2" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar05_Tau2', 'CONEPASYM', -100.)".format(stream, line),
			"CONEPTASYM_tau2" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar05_Tau2', 'CONEPTASYM', -100.)".format(stream, line),
			"CONEDELTAETA_tau2" : "RELINFO('/Event/{0}/Phys/{1}/P2ConeVar05_Tau2', 'CONEDELTAETA', -100.)".format(stream, line)
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
def addKTTDTF(branch, int):
    dtf = branch.addTupleTool("B2KtautauDTF1/ConsBp_{}".format(int))
    dtf.daughtersToConstrain = ['nu_tau', 'nu_tau~', 'tau-', 'tau+']
    dtf.UseFullTreeInName = True
    dtf.UpdateDaughters = True
    dtf.Verbose = True
    dtf.constrainToOriginVertex = True
    dtf.initStrategy = int

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

def addTools( DecayTreeTuple, stream, line, year ):
    #run modified DTF
    addKTTDTF(DecayTreeTuple.Bp, 0) # default
    addKTTDTF(DecayTreeTuple.Bp, 1) # pions direction
    #add some helpful variables for B0
    addTopInfo(DecayTreeTuple.Bp)
    #add isolation info
    addIsoInfo(DecayTreeTuple, stream, line, year)
    #add submasses
    addSubmassInfo(DecayTreeTuple.Bp)
    #add helpful info for final state tracks
    addTrkInfo(DecayTreeTuple.taup_pi1)
    addTrkInfo(DecayTreeTuple.taup_pi2)
    addTrkInfo(DecayTreeTuple.taup_pi3)
    addTrkInfo(DecayTreeTuple.taum_pi1)
    addTrkInfo(DecayTreeTuple.taum_pi2)
    addTrkInfo(DecayTreeTuple.taum_pi3)
    addTrkInfo(DecayTreeTuple.Kp)
    #add M12, M23, M13 for tau's
    addMInfo(DecayTreeTuple.taup)
    addMInfo(DecayTreeTuple.taum)
    # add trigger information
    addTISTOS(DecayTreeTuple.Bp)

# Stream and stripping line we want to use
if isMC:
	stream = 'AllStreams'
else:	
	stream = 'Bhadron'
line = 'B2KTauTauLine'
line_SS = 'B2KTauTauSSLine'

# Create an ntuple to capture B+ decays from the StrippingLine line
dtt = DecayTreeTuple('ntuple')
dtt_SS = DecayTreeTuple('ntuple_SS')
mcdtt = MCDecayTreeTuple('mc_ntuple')

#Loose preselections to reduce the size of tuples
if isMC:
	b_strip = AutomaticData('/Event/{0}/Phys/{1}/Particles'.format(stream, line))# DST
else:
	b_strip    = AutomaticData('/Phys/{0}/Particles'.format(line))# microDST 
	b_strip_SS = AutomaticData('/Phys/{0}/Particles'.format(line_SS))# microDST 

filtercode = " (NINGENERATION( (M < 3400) & (ABSID=='B+') , 0) == 0)  "\
			 "& (NINGENERATION( (M > 5000) & (ABSID=='B+') , 0) == 0) "\
			 "& (NINGENERATION( (M < 800) & (ABSID=='tau+') , 1) == 0) "\
			 "& (NINGENERATION( (M > 1600) & (ABSID=='tau+') , 1) == 0) "

b_Filter    = FilterSelection("b_Filter", b_strip, Code = filtercode)
if not isMC:
	b_Filter_SS = FilterSelection("b_Filter_SS", b_strip_SS, Code = filtercode)

if isMC:
	dtt.Inputs = [b_Filter.outputLocation()]
	#dtt.Inputs = ['/Event/{0}/Phys/{1}/Particles'.format(stream, line)] # DST
else:
	dtt.Inputs = [b_Filter.outputLocation()]
	dtt_SS.Inputs = [b_Filter_SS.outputLocation()]
	# dtt.Inputs = ['/Phys/{0}/Particles'.format(line)] # microDST 
	# dtt_SS.Inputs = ['/Phys/{0}/Particles'.format(line_SS)] # microDST 

dtt.setDescriptorTemplate('${Bp}[B+ -> ${taup}(tau+ -> ${taup_pi1}pi+ ${taup_pi2}pi- ${taup_pi3}pi+) ${taum}(tau- -> ${taum_pi1}pi- ${taum_pi2}pi+ ${taum_pi3}pi-) ${Kp}K+]CC')
dtt_SS.setDescriptorTemplate('(${Bp}[B+ -> ${taup}(tau+ -> ${taup_pi1}pi+ ${taup_pi2}pi- ${taup_pi3}pi+) ${taum}(tau+ -> ${taum_pi1}pi+ ${taum_pi2}pi- ${taum_pi3}pi+) ${Kp}K+]CC)'
						  '|| (${Bp}[B+ -> ${taup}(tau+ -> ${taup_pi1}pi+ ${taup_pi2}pi- ${taup_pi3}pi+) ${taum}(tau+ -> ${taum_pi1}pi+ ${taum_pi2}pi- ${taum_pi3}pi+) ${Kp}K-]CC)')
# Generated (Acceptance cut: DaughtersInLHCb)
mcdtt.Decay = '[B+ => ^(D~0 => ^(a_1(1260)- => ^(rho(770)0 => ^pi+ ^pi-) ^pi-) ^K+)  ^(D+ => ^(a_1(1260)+ => ^(rho(770)0 => ^pi+ ^pi-) ^pi+) ^K0) ]CC' # add {X0} for tau decay

# TUPLE TOOLS
ToolList = ['TupleToolTrackInfo', 'TupleToolVtxIsoln', 'TupleToolPrimaries', 'TupleToolCorrectedMass', 'TupleToolPropertime', 'TupleToolAngles', 'TupleToolL0Calo', 'TupleToolRecoStats']
dtt.ToolList += ToolList
if not isMC:
	dtt_SS.ToolList += ToolList
if isMC:
	dtt.addTupleTool('TupleToolMCTruth')
	dtt.addTupleTool('TupleToolMCBackgroundInfo')

# Adds several tupletools, including DTF
addTools(dtt, stream, line, year)
if not isMC:
	addTools(dtt_SS, stream, line_SS, year)

GENToolList = ['MCTupleToolKinematic', 'MCTupleToolHierarchy', 'TupleToolEventInfo']
if isMC:
	mcdtt.ToolList += GENToolList

from PhysConf.Filters import LoKi_Filters
fltrs = LoKi_Filters (
    STRIP_Code = "HLT_PASS_RE('StrippingB2KTauTauLineDecision')"
)
if not isMC:
	DaVinci().EventPreFilters = fltrs.filters('Filters')

# PV checker -> TupleToolPrimaries requires there to be a PV in the event, otherwise it throws an error
DaVinci().UserAlgorithms += [CheckPV()]

DaVinci().UserAlgorithms += [b_Filter]
if not isMC:
	DaVinci().UserAlgorithms += [b_Filter_SS]

DaVinci().UserAlgorithms += [dtt]
if isMC:
	DaVinci().UserAlgorithms += [mcdtt]
if not isMC:
	DaVinci().UserAlgorithms += [dtt_SS]

# Configure DaVinci
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

IOHelper().inputFiles(['/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135702.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135705.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135706.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135707.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135709.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135710.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135712.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135715.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135717.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135718.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135719.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135720.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135721.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135722.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135723.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135725.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135726.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135727.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135728.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135730.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135733.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135736.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135737.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135740.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135741.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135742.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135743.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135745.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135746.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135747.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135748.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_902135750.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095902.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095903.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095904.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095905.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095906.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095907.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095908.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095909.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095910.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095911.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095912.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095914.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095915.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095917.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095918.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095919.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095922.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095924.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095925.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095927.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095928.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095929.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095931.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095932.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095933.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095934.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095935.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095936.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095937.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095938.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095939.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095941.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095942.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095943.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095947.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095949.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/1000evts_s28r2_906095950.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095104.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095112.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095128.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095143.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095152.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095165.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095195.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095196.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095200.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095208.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095209.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095217.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095243.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095245.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095248.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095250.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095262.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095267.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095288.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095289.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095293.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095303.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095320.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095350.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095355.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095460.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095540.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095558.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095570.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095607.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095649.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902095653.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902135725.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902135746.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902135748.dst',
					'/panfs/felician/SimulationJobs/12197098/2016/Sim09h_ReDecay/MagUp/100evts_s28r2_902135750.dst'])
