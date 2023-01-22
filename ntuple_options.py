from Configurables import DecayTreeTuple
from DecayTreeTuple.Configuration import *
from Configurables import LoKi__Hybrid__TupleTool
from Configurables import CheckPV
from Configurables import DaVinci
from PhysConf.Selections import FilterSelection
from PhysSelPython.Wrappers import AutomaticData

isMC = True
year = '2016'
pol = 'MagDown'

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
def addKTTDTF(branch):
	# both taus directions are initialised based on vertices
    dtf0 = branch.addTupleTool("B2KtautauDTF1/ConsBp_0")
    dtf0.daughtersToConstrain = ['nu_tau', 'nu_tau~', 'tau-', 'tau+']
    dtf0.UseFullTreeInName = True
    dtf0.UpdateDaughters = True
    dtf0.Verbose = True
    dtf0.constrainToOriginVertex = True
    dtf0.initStrategy = 0

	# both taus directions are initialised based on visible 3pi momenta
    dtf1 = branch.addTupleTool("B2KtautauDTF1/ConsBp_1")
    dtf1.daughtersToConstrain = ['nu_tau', 'nu_tau~', 'tau-', 'tau+']
    dtf1.UseFullTreeInName = True
    dtf1.UpdateDaughters = True
    dtf1.Verbose = True
    dtf1.constrainToOriginVertex = True
    dtf1.initStrategy = 1

	# sequential DTF (strategy0 -> strategy1 -> strategy2: using Marseille's tau momentum)
    dtf_seq = branch.addTupleTool("TupleToolKTauTauDTFSequential/ConsBp_seq")
    dtf_seq.daughtersToConstrain = ['nu_tau', 'nu_tau~', 'tau-', 'tau+']
    dtf_seq.UseFullTreeInName = True
    dtf_seq.UpdateDaughters = True
    dtf_seq.Verbose = True
    dtf_seq.constrainToOriginVertex = True


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
    addKTTDTF(DecayTreeTuple.Bp) # default
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
mcdtt.Decay = '[B+ => ^(tau+ ==> ^pi+ ^pi- ^pi+ ^nu_tau~) ^(tau- ==> ^pi- ^pi+ ^pi- ^nu_tau) K+]CC' # add {X0} for tau decay

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
	DaVinci().TupleFile = '/panfs/felician/B2Ktautau/ROOT_Sim/'+year+'/DVntuple_MC_'+year+'_'+pol+'.root'
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
		if year == '2016':
			DaVinci().CondDBtag = "sim-20170721-2-vc-mu100"
		elif year == '2017':
			DaVinci().CondDBtag = "sim-20190430-1-vc-mu100" 
		elif year == '2018':       
			DaVinci().CondDBtag = "sim-20190430-vc-mu100"   
	else:
		if year == '2016':
			DaVinci().CondDBtag = "sim-20170721-2-vc-md100"
		elif year == '2017':
			DaVinci().CondDBtag = "sim-20190430-1-vc-md100"
		elif year == '2018':
			DaVinci().CondDBtag = "sim-20190430-vc-md100"

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

IOHelper().inputFiles(['/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143251.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143252.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143253.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143256.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143257.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143258.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143260.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143264.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143265.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143266.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143267.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143268.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143269.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143270.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143271.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143272.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143275.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143276.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143277.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143278.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143279.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143280.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143282.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143283.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143284.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143285.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143287.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143289.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143290.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143291.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143292.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143293.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1011143299.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144255.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144258.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144259.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144261.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144262.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144264.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144265.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144266.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144268.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144269.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144270.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144271.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144273.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144274.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144275.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144276.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144278.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144279.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144280.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144281.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144282.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144283.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144285.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144286.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144287.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144288.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144289.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144290.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144292.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144294.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144295.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144297.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144298.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1013144300.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1016090662.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1016090663.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1016090664.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1016090665.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1016090668.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1016090669.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1016090670.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1016090674.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1016090675.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1016090681.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1016090687.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1016090696.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1016090700.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085751.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085752.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085753.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085754.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085755.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085760.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085761.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085762.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085763.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085764.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085766.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085767.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085768.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085769.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085770.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085773.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085774.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085775.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085776.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085777.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085778.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085779.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085780.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085781.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085782.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085785.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085787.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085788.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085789.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085790.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085793.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085794.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085796.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085797.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085800.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085801.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085802.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085803.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085807.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085808.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085809.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085810.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085812.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085813.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085814.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085815.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085818.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085819.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085820.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085827.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085831.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085832.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085834.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085837.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085838.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085839.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085840.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085841.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085842.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085846.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085847.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085848.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085849.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085850.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085851.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085853.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085855.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085858.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085859.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085862.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085864.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085866.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085868.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085869.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085871.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085872.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085873.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085875.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085877.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085878.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085880.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085881.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085882.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085884.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085885.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085886.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085887.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085888.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085889.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085890.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085891.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085892.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085893.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085894.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085895.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085897.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085899.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085900.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085902.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085905.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085906.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085907.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085908.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085909.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085910.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085912.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085914.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085915.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085919.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085921.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085922.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085923.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085927.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085929.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085930.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085931.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085934.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085935.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085936.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085937.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085938.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085939.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085940.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085941.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085942.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085943.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085944.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085945.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085946.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085948.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085950.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085951.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085953.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085954.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085955.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085956.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085957.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085958.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085959.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085960.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085961.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085962.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085963.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085965.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085966.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085967.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085968.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085969.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085970.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085971.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085972.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085973.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085975.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085977.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085979.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085980.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085982.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085984.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085986.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085987.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085988.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085990.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085991.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085992.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085994.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085995.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019085998.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_1019086000.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102551.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102552.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102554.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102555.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102556.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102557.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102559.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102560.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102561.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102562.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102565.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102568.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102569.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102570.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102571.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102573.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102574.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102575.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102576.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102579.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102580.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102581.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102582.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102584.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102585.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102586.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102589.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102590.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102591.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102592.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102593.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102594.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102596.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102597.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102598.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/1000evts_s28r2_908102600.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_1007162910.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_1007162926.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_1007162935.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_1007162947.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_1007162958.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_1007162981.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_1007163003.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_1007163022.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_1007163074.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_1007163085.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_1007163121.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_1007163144.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_1007163170.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_1007163180.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_1007163201.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_1007163205.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_1007163209.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_1007163222.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_1007163223.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_1007163231.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_1007163283.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_1007163289.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_1007163322.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_1007163338.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_1007163341.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_1007163373.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_1007163377.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_1007163386.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/100evts_s28r2_1007163388.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/2000evts_s28r2_1018164952.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/2000evts_s28r2_1018164955.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/2000evts_s28r2_1018164958.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/2000evts_s28r2_1018164960.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/2000evts_s28r2_1018164961.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/2000evts_s28r2_1018164962.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/2000evts_s28r2_1018164963.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/2000evts_s28r2_1018164964.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/2000evts_s28r2_1018164965.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/2000evts_s28r2_1018164966.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/2000evts_s28r2_1018164967.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/2000evts_s28r2_1018164968.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/2000evts_s28r2_1018164969.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/2000evts_s28r2_1018164970.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/2000evts_s28r2_1018164980.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/2000evts_s28r2_1018164981.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/2000evts_s28r2_1018164982.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/2000evts_s28r2_1018164983.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/2000evts_s28r2_1018164984.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/2000evts_s28r2_1018164986.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/2000evts_s28r2_1018164987.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/2000evts_s28r2_1018164988.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/2000evts_s28r2_1018164989.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/2000evts_s28r2_1018164991.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/2000evts_s28r2_1018164993.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/2000evts_s28r2_1018164994.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/2000evts_s28r2_1018164995.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/2000evts_s28r2_1018164996.dst',
					'/panfs/felician/SimulationJobs/12101010/2016/Sim09h_ReDecay/MagDown/2000evts_s28r2_1018164998.dst'
 					], clear=True)
