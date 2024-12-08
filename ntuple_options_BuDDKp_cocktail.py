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
applySel = False
if(not isMC):
	applySel = True

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
		"M13" : "M13",
		"AMAXDOCA": "PFUNA(AMAXDOCA('LoKi::TrgDistanceCalculator'))",
		"AMINDOCA": "PFUNA(AMINDOCA('LoKi::TrgDistanceCalculator'))",
		"DOCACHI2MAX" : "DOCACHI2MAX"
}

def addSubmassInfo( branch ):
	SubmassInfo = branch.addTupleTool("TupleToolSubMass")
	SubmassInfo.SetMax = 3
	SubmassInfo.SubVertexFit = False

# Isolation information (stripping line output isolation variables)
def addIsoInfo(branch, stream, line, year):
	LoKi_Cone =  branch.Bp.addTupleTool("LoKi::Hybrid::TupleTool/LoKi_Cone")

	LoKi_Cone.Variables = {
		# Vertex isolation variables
		# B+
		"VTXISONUMVTX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_VertexIsoInfo', 'VTXISONUMVTX', -100)".format(stream,line),
		"VTXISODCHI2ONETRACK_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_VertexIsoInfo', 'VTXISODCHI2ONETRACK', -100)".format(stream,line),
		"VTXISODCHI2MASSONETRACK_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_VertexIsoInfo', 'VTXISODCHI2MASSONETRACK', -100)".format(stream,line),
		"VTXISODCHI2TWOTRACK_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_VertexIsoInfo', 'VTXISODCHI2TWOTRACK', -100)".format(stream,line),
		"VTXISODCHI2MASSTWOTRACK_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_VertexIsoInfo', 'VTXISODCHI2MASSTWOTRACK', -100)".format(stream,line),

		# tau+
		"VTXISONUMVTX_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_VertexIsoInfo', 'VTXISONUMVTX', -100)".format(stream,line),
		"VTXISODCHI2ONETRACK_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_VertexIsoInfo', 'VTXISODCHI2ONETRACK', -100)".format(stream,line),
		"VTXISODCHI2MASSONETRACK_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_VertexIsoInfo', 'VTXISODCHI2MASSONETRACK', -100)".format(stream,line),
		"VTXISODCHI2TWOTRACK_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_VertexIsoInfo', 'VTXISODCHI2TWOTRACK', -100)".format(stream,line),
		"VTXISODCHI2MASSTWOTRACK_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_VertexIsoInfo', 'VTXISODCHI2MASSTWOTRACK', -100)".format(stream,line),

		# tau-
		"VTXISONUMVTX_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_VertexIsoInfo', 'VTXISONUMVTX', -100)".format(stream,line),
		"VTXISODCHI2ONETRACK_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_VertexIsoInfo', 'VTXISODCHI2ONETRACK', -100)".format(stream,line),
		"VTXISODCHI2MASSONETRACK_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_VertexIsoInfo', 'VTXISODCHI2MASSONETRACK', -100)".format(stream,line),
		"VTXISODCHI2TWOTRACK_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_VertexIsoInfo', 'VTXISODCHI2TWOTRACK', -100)".format(stream,line),
		"VTXISODCHI2MASSTWOTRACK_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_VertexIsoInfo', 'VTXISODCHI2MASSTWOTRACK', -100)".format(stream,line),

		# Come isolation variables
		# B+ (cone size = 0.5)
		"CC_05_ANGLE_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'CC_ANGLE', -100.)".format(stream, line),
		"CC_05_MULT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'CC_MULT', -100.)".format(stream, line),
		"CC_05_SPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'CC_SPT', -100.)".format(stream, line),
		"CC_05_VPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'CC_VPT', -100.)".format(stream, line),
		"CC_05_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'CC_PX', -100.)".format(stream, line),
		"CC_05_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'CC_PY', -100.)".format(stream, line),
		"CC_05_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'CC_PZ', -100.)".format(stream, line),
		"CC_05_PASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'CC_PASYM', -100.)".format(stream, line),
		"CC_05_PTASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'CC_PTASYM', -100.)".format(stream, line),
		"CC_05_PXASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'CC_PXASYM', -100.)".format(stream, line),
		"CC_05_PYASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'CC_PYASYM', -100.)".format(stream, line),
		"CC_05_PZASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'CC_PZASYM', -100.)".format(stream, line),
		"CC_05_DELTAETA_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'CC_DELTAETA', -100.)".format(stream, line),
		"CC_05_DELTAPHI_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'CC_DELTAPHI', -100.)".format(stream, line),
		"CC_05_IT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'CC_IT', -100.)".format(stream, line),
		"CC_05_MAXPT_Q_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'CC_MAXPT_Q', -100.)".format(stream, line),
		"CC_05_MAXPT_PT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'CC_MAXPT_PT', -100.)".format(stream, line),
		"CC_05_MAXPT_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'CC_MAXPT_PX', -100.)".format(stream, line),
		"CC_05_MAXPT_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'CC_MAXPT_PY', -100.)".format(stream, line),
		"CC_05_MAXPT_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'CC_MAXPT_PZ', -100.)".format(stream, line),
		"CC_05_MAXPT_PE_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'CC_MAXPT_PE', -100.)".format(stream, line),
		"NC_05_ANGLE_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'NC_ANGLE', -100.)".format(stream, line),
		"NC_05_MULT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'NC_MULT', -100.)".format(stream, line),
		"NC_05_SPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'NC_SPT', -100.)".format(stream, line),
		"NC_05_VPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'NC_VPT', -100.)".format(stream, line),
		"NC_05_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'NC_PX', -100.)".format(stream, line),
		"NC_05_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'NC_PY', -100.)".format(stream, line),
		"NC_05_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'NC_PZ', -100.)".format(stream, line),
		"NC_05_PASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'NC_PASYM', -100.)".format(stream, line),
		"NC_05_PTASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'NC_PTASYM', -100.)".format(stream, line),
		"NC_05_PXASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'NC_PXASYM', -100.)".format(stream, line),
		"NC_05_PYASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'NC_PYASYM', -100.)".format(stream, line),
		"NC_05_PZASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'NC_PZASYM', -100.)".format(stream, line),
		"NC_05_DELTAETA_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'NC_DELTAETA', -100.)".format(stream, line),
		"NC_05_DELTAPHI_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'NC_DELTAPHI', -100.)".format(stream, line),
		"NC_05_IT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'NC_IT', -100.)".format(stream, line),
		"NC_05_MAXPT_PT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'NC_MAXPT_PT', -100.)".format(stream, line),
		"NC_05_MAXPT_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'NC_MAXPT_PX', -100.)".format(stream, line),
		"NC_05_MAXPT_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'NC_MAXPT_PY', -100.)".format(stream, line),
		"NC_05_MAXPT_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'NC_MAXPT_PZ', -100.)".format(stream, line),
		"CCNC_05_IT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'CCNC_IT', -100.)".format(stream, line),

		# K+ (cone size = 0.5)
		"CC_05_ANGLE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'CC_ANGLE', -100.)".format(stream, line),
		"CC_05_MULT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'CC_MULT', -100.)".format(stream, line),
		"CC_05_SPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'CC_SPT', -100.)".format(stream, line),
		"CC_05_VPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'CC_VPT', -100.)".format(stream, line),
		"CC_05_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'CC_PX', -100.)".format(stream, line),
		"CC_05_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'CC_PY', -100.)".format(stream, line),
		"CC_05_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'CC_PZ', -100.)".format(stream, line),
		"CC_05_PASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'CC_PASYM', -100.)".format(stream, line),
		"CC_05_PTASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'CC_PTASYM', -100.)".format(stream, line),
		"CC_05_PXASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'CC_PXASYM', -100.)".format(stream, line),
		"CC_05_PYASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'CC_PYASYM', -100.)".format(stream, line),
		"CC_05_PZASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'CC_PZASYM', -100.)".format(stream, line),
		"CC_05_DELTAETA_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'CC_DELTAETA', -100.)".format(stream, line),
		"CC_05_DELTAPHI_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'CC_DELTAPHI', -100.)".format(stream, line),
		"CC_05_IT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'CC_IT', -100.)".format(stream, line),
		"CC_05_MAXPT_Q_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'CC_MAXPT_Q', -100.)".format(stream, line),
		"CC_05_MAXPT_PT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'CC_MAXPT_PT', -100.)".format(stream, line),
		"CC_05_MAXPT_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'CC_MAXPT_PX', -100.)".format(stream, line),
		"CC_05_MAXPT_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'CC_MAXPT_PY', -100.)".format(stream, line),
		"CC_05_MAXPT_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'CC_MAXPT_PZ', -100.)".format(stream, line),
		"CC_05_MAXPT_PE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'CC_MAXPT_PE', -100.)".format(stream, line),
		"NC_05_ANGLE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'NC_ANGLE', -100.)".format(stream, line),
		"NC_05_MULT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'NC_MULT', -100.)".format(stream, line),
		"NC_05_SPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'NC_SPT', -100.)".format(stream, line),
		"NC_05_VPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'NC_VPT', -100.)".format(stream, line),
		"NC_05_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'NC_PX', -100.)".format(stream, line),
		"NC_05_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'NC_PY', -100.)".format(stream, line),
		"NC_05_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'NC_PZ', -100.)".format(stream, line),
		"NC_05_PASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'NC_PASYM', -100.)".format(stream, line),
		"NC_05_PTASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'NC_PTASYM', -100.)".format(stream, line),
		"NC_05_PXASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'NC_PXASYM', -100.)".format(stream, line),
		"NC_05_PYASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'NC_PYASYM', -100.)".format(stream, line),
		"NC_05_PZASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'NC_PZASYM', -100.)".format(stream, line),
		"NC_05_DELTAETA_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'NC_DELTAETA', -100.)".format(stream, line),
		"NC_05_DELTAPHI_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'NC_DELTAPHI', -100.)".format(stream, line),
		"NC_05_IT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'NC_IT', -100.)".format(stream, line),
		"NC_05_MAXPT_PT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'NC_MAXPT_PT', -100.)".format(stream, line),
		"NC_05_MAXPT_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'NC_MAXPT_PX', -100.)".format(stream, line),
		"NC_05_MAXPT_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'NC_MAXPT_PY', -100.)".format(stream, line),
		"NC_05_MAXPT_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'NC_MAXPT_PZ', -100.)".format(stream, line),
		"CCNC_05_IT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'CCNC_IT', -100.)".format(stream, line),

		# tau+ (cone size = 0.5)
		"CC_05_ANGLE_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'CC_ANGLE', -100.)".format(stream, line),
		"CC_05_MULT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'CC_MULT', -100.)".format(stream, line),
		"CC_05_SPT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'CC_SPT', -100.)".format(stream, line),
		"CC_05_VPT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'CC_VPT', -100.)".format(stream, line),
		"CC_05_PX_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'CC_PX', -100.)".format(stream, line),
		"CC_05_PY_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'CC_PY', -100.)".format(stream, line),
		"CC_05_PZ_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'CC_PZ', -100.)".format(stream, line),
		"CC_05_PASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'CC_PASYM', -100.)".format(stream, line),
		"CC_05_PTASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'CC_PTASYM', -100.)".format(stream, line),
		"CC_05_PXASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'CC_PXASYM', -100.)".format(stream, line),
		"CC_05_PYASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'CC_PYASYM', -100.)".format(stream, line),
		"CC_05_PZASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'CC_PZASYM', -100.)".format(stream, line),
		"CC_05_DELTAETA_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'CC_DELTAETA', -100.)".format(stream, line),
		"CC_05_DELTAPHI_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'CC_DELTAPHI', -100.)".format(stream, line),
		"CC_05_IT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'CC_IT', -100.)".format(stream, line),
		"CC_05_MAXPT_Q_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'CC_MAXPT_Q', -100.)".format(stream, line),
		"CC_05_MAXPT_PT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'CC_MAXPT_PT', -100.)".format(stream, line),
		"CC_05_MAXPT_PX_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'CC_MAXPT_PX', -100.)".format(stream, line),
		"CC_05_MAXPT_PY_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'CC_MAXPT_PY', -100.)".format(stream, line),
		"CC_05_MAXPT_PZ_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'CC_MAXPT_PZ', -100.)".format(stream, line),
		"CC_05_MAXPT_PE_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'CC_MAXPT_PE', -100.)".format(stream, line),
		"NC_05_ANGLE_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'NC_ANGLE', -100.)".format(stream, line),
		"NC_05_MULT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'NC_MULT', -100.)".format(stream, line),
		"NC_05_SPT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'NC_SPT', -100.)".format(stream, line),
		"NC_05_VPT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'NC_VPT', -100.)".format(stream, line),
		"NC_05_PX_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'NC_PX', -100.)".format(stream, line),
		"NC_05_PY_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'NC_PY', -100.)".format(stream, line),
		"NC_05_PZ_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'NC_PZ', -100.)".format(stream, line),
		"NC_05_PASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'NC_PASYM', -100.)".format(stream, line),
		"NC_05_PTASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'NC_PTASYM', -100.)".format(stream, line),
		"NC_05_PXASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'NC_PXASYM', -100.)".format(stream, line),
		"NC_05_PYASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'NC_PYASYM', -100.)".format(stream, line),
		"NC_05_PZASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'NC_PZASYM', -100.)".format(stream, line),
		"NC_05_DELTAETA_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'NC_DELTAETA', -100.)".format(stream, line),
		"NC_05_DELTAPHI_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'NC_DELTAPHI', -100.)".format(stream, line),
		"NC_05_IT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'NC_IT', -100.)".format(stream, line),
		"NC_05_MAXPT_PT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'NC_MAXPT_PT', -100.)".format(stream, line),
		"NC_05_MAXPT_PX_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'NC_MAXPT_PX', -100.)".format(stream, line),
		"NC_05_MAXPT_PY_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'NC_MAXPT_PY', -100.)".format(stream, line),
		"NC_05_MAXPT_PZ_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'NC_MAXPT_PZ', -100.)".format(stream, line),
		"CCNC_05_IT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone05', 'CCNC_IT', -100.)".format(stream, line),

		# tau- (cone size = 0.5)
		"CC_05_ANGLE_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'CC_ANGLE', -100.)".format(stream, line),
		"CC_05_MULT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'CC_MULT', -100.)".format(stream, line),
		"CC_05_SPT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'CC_SPT', -100.)".format(stream, line),
		"CC_05_VPT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'CC_VPT', -100.)".format(stream, line),
		"CC_05_PX_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'CC_PX', -100.)".format(stream, line),
		"CC_05_PY_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'CC_PY', -100.)".format(stream, line),
		"CC_05_PZ_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'CC_PZ', -100.)".format(stream, line),
		"CC_05_PASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'CC_PASYM', -100.)".format(stream, line),
		"CC_05_PTASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'CC_PTASYM', -100.)".format(stream, line),
		"CC_05_PXASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'CC_PXASYM', -100.)".format(stream, line),
		"CC_05_PYASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'CC_PYASYM', -100.)".format(stream, line),
		"CC_05_PZASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'CC_PZASYM', -100.)".format(stream, line),
		"CC_05_DELTAETA_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'CC_DELTAETA', -100.)".format(stream, line),
		"CC_05_DELTAPHI_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'CC_DELTAPHI', -100.)".format(stream, line),
		"CC_05_IT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'CC_IT', -100.)".format(stream, line),
		"CC_05_MAXPT_Q_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'CC_MAXPT_Q', -100.)".format(stream, line),
		"CC_05_MAXPT_PT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'CC_MAXPT_PT', -100.)".format(stream, line),
		"CC_05_MAXPT_PX_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'CC_MAXPT_PX', -100.)".format(stream, line),
		"CC_05_MAXPT_PY_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'CC_MAXPT_PY', -100.)".format(stream, line),
		"CC_05_MAXPT_PZ_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'CC_MAXPT_PZ', -100.)".format(stream, line),
		"CC_05_MAXPT_PE_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'CC_MAXPT_PE', -100.)".format(stream, line),
		"NC_05_ANGLE_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'NC_ANGLE', -100.)".format(stream, line),
		"NC_05_MULT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'NC_MULT', -100.)".format(stream, line),
		"NC_05_SPT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'NC_SPT', -100.)".format(stream, line),
		"NC_05_VPT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'NC_VPT', -100.)".format(stream, line),
		"NC_05_PX_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'NC_PX', -100.)".format(stream, line),
		"NC_05_PY_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'NC_PY', -100.)".format(stream, line),
		"NC_05_PZ_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'NC_PZ', -100.)".format(stream, line),
		"NC_05_PASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'NC_PASYM', -100.)".format(stream, line),
		"NC_05_PTASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'NC_PTASYM', -100.)".format(stream, line),
		"NC_05_PXASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'NC_PXASYM', -100.)".format(stream, line),
		"NC_05_PYASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'NC_PYASYM', -100.)".format(stream, line),
		"NC_05_PZASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'NC_PZASYM', -100.)".format(stream, line),
		"NC_05_DELTAETA_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'NC_DELTAETA', -100.)".format(stream, line),
		"NC_05_DELTAPHI_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'NC_DELTAPHI', -100.)".format(stream, line),
		"NC_05_IT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'NC_IT', -100.)".format(stream, line),
		"NC_05_MAXPT_PT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'NC_MAXPT_PT', -100.)".format(stream, line),
		"NC_05_MAXPT_PX_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'NC_MAXPT_PX', -100.)".format(stream, line),
		"NC_05_MAXPT_PY_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'NC_MAXPT_PY', -100.)".format(stream, line),
		"NC_05_MAXPT_PZ_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'NC_MAXPT_PZ', -100.)".format(stream, line),
		"CCNC_05_IT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone05', 'CCNC_IT', -100.)".format(stream, line),

		# B+ (cone size = 1.0)
		# "CC_10_ANGLE_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_ANGLE', -100.)".format(stream, line),
		# "CC_10_MULT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_MULT', -100.)".format(stream, line),
		# "CC_10_SPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_SPT', -100.)".format(stream, line),
		# "CC_10_VPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_VPT', -100.)".format(stream, line),
		# "CC_10_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_PX', -100.)".format(stream, line),
		# "CC_10_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_PY', -100.)".format(stream, line),
		# "CC_10_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_PZ', -100.)".format(stream, line),
		# "CC_10_PASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_PASYM', -100.)".format(stream, line),
		# "CC_10_PTASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_PTASYM', -100.)".format(stream, line),
		# "CC_10_PXASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_PXASYM', -100.)".format(stream, line),
		# "CC_10_PYASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_PYASYM', -100.)".format(stream, line),
		# "CC_10_PZASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_PZASYM', -100.)".format(stream, line),
		# "CC_10_DELTAETA_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_DELTAETA', -100.)".format(stream, line),
		# "CC_10_DELTAPHI_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_DELTAPHI', -100.)".format(stream, line),
		# "CC_10_IT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_IT', -100.)".format(stream, line),
		# "CC_10_MAXPT_Q_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_MAXPT_Q', -100.)".format(stream, line),
		# "CC_10_MAXPT_PT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_MAXPT_PT', -100.)".format(stream, line),
		# "CC_10_MAXPT_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_MAXPT_PX', -100.)".format(stream, line),
		# "CC_10_MAXPT_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_MAXPT_PY', -100.)".format(stream, line),
		# "CC_10_MAXPT_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_MAXPT_PZ', -100.)".format(stream, line),
		# "CC_10_MAXPT_PE_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_MAXPT_PE', -100.)".format(stream, line),
		# "NC_10_ANGLE_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_ANGLE', -100.)".format(stream, line),
		# "NC_10_MULT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_MULT', -100.)".format(stream, line),
		# "NC_10_SPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_SPT', -100.)".format(stream, line),
		# "NC_10_VPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_VPT', -100.)".format(stream, line),
		# "NC_10_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_PX', -100.)".format(stream, line),
		# "NC_10_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_PY', -100.)".format(stream, line),
		# "NC_10_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_PZ', -100.)".format(stream, line),
		# "NC_10_PASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_PASYM', -100.)".format(stream, line),
		# "NC_10_PTASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_PTASYM', -100.)".format(stream, line),
		# "NC_10_PXASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_PXASYM', -100.)".format(stream, line),
		# "NC_10_PYASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_PYASYM', -100.)".format(stream, line),
		# "NC_10_PZASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_PZASYM', -100.)".format(stream, line),
		# "NC_10_DELTAETA_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_DELTAETA', -100.)".format(stream, line),
		# "NC_10_DELTAPHI_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_DELTAPHI', -100.)".format(stream, line),
		# "NC_10_IT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_IT', -100.)".format(stream, line),
		# "NC_10_MAXPT_PT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_MAXPT_PT', -100.)".format(stream, line),
		# "NC_10_MAXPT_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_MAXPT_PX', -100.)".format(stream, line),
		# "NC_10_MAXPT_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_MAXPT_PY', -100.)".format(stream, line),
		# "NC_10_MAXPT_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_MAXPT_PZ', -100.)".format(stream, line),
		# "CCNC_10_IT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CCNC_IT', -100.)".format(stream, line),

		# # K+ (cone size = 1.0)
		# "CC_10_ANGLE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_ANGLE', -100.)".format(stream, line),
		# "CC_10_MULT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_MULT', -100.)".format(stream, line),
		# "CC_10_SPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_SPT', -100.)".format(stream, line),
		# "CC_10_VPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_VPT', -100.)".format(stream, line),
		# "CC_10_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_PX', -100.)".format(stream, line),
		# "CC_10_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_PY', -100.)".format(stream, line),
		# "CC_10_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_PZ', -100.)".format(stream, line),
		# "CC_10_PASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_PASYM', -100.)".format(stream, line),
		# "CC_10_PTASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_PTASYM', -100.)".format(stream, line),
		# "CC_10_PXASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_PXASYM', -100.)".format(stream, line),
		# "CC_10_PYASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_PYASYM', -100.)".format(stream, line),
		# "CC_10_PZASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_PZASYM', -100.)".format(stream, line),
		# "CC_10_DELTAETA_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_DELTAETA', -100.)".format(stream, line),
		# "CC_10_DELTAPHI_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_DELTAPHI', -100.)".format(stream, line),
		# "CC_10_IT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_IT', -100.)".format(stream, line),
		# "CC_10_MAXPT_Q_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_MAXPT_Q', -100.)".format(stream, line),
		# "CC_10_MAXPT_PT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_MAXPT_PT', -100.)".format(stream, line),
		# "CC_10_MAXPT_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_MAXPT_PX', -100.)".format(stream, line),
		# "CC_10_MAXPT_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_MAXPT_PY', -100.)".format(stream, line),
		# "CC_10_MAXPT_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_MAXPT_PZ', -100.)".format(stream, line),
		# "CC_10_MAXPT_PE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_MAXPT_PE', -100.)".format(stream, line),
		# "NC_10_ANGLE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_ANGLE', -100.)".format(stream, line),
		# "NC_10_MULT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_MULT', -100.)".format(stream, line),
		# "NC_10_SPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_SPT', -100.)".format(stream, line),
		# "NC_10_VPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_VPT', -100.)".format(stream, line),
		# "NC_10_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_PX', -100.)".format(stream, line),
		# "NC_10_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_PY', -100.)".format(stream, line),
		# "NC_10_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_PZ', -100.)".format(stream, line),
		# "NC_10_PASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_PASYM', -100.)".format(stream, line),
		# "NC_10_PTASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_PTASYM', -100.)".format(stream, line),
		# "NC_10_PXASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_PXASYM', -100.)".format(stream, line),
		# "NC_10_PYASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_PYASYM', -100.)".format(stream, line),
		# "NC_10_PZASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_PZASYM', -100.)".format(stream, line),
		# "NC_10_DELTAETA_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_DELTAETA', -100.)".format(stream, line),
		# "NC_10_DELTAPHI_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_DELTAPHI', -100.)".format(stream, line),
		# "NC_10_IT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_IT', -100.)".format(stream, line),
		# "NC_10_MAXPT_PT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_MAXPT_PT', -100.)".format(stream, line),
		# "NC_10_MAXPT_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_MAXPT_PX', -100.)".format(stream, line),
		# "NC_10_MAXPT_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_MAXPT_PY', -100.)".format(stream, line),
		# "NC_10_MAXPT_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_MAXPT_PZ', -100.)".format(stream, line),
		# "CCNC_10_IT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CCNC_IT', -100.)".format(stream, line),

		# # tau+ (cone size = 1.0)
		# "CC_10_ANGLE_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'CC_ANGLE', -100.)".format(stream, line),
		# "CC_10_MULT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'CC_MULT', -100.)".format(stream, line),
		# "CC_10_SPT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'CC_SPT', -100.)".format(stream, line),
		# "CC_10_VPT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'CC_VPT', -100.)".format(stream, line),
		# "CC_10_PX_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'CC_PX', -100.)".format(stream, line),
		# "CC_10_PY_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'CC_PY', -100.)".format(stream, line),
		# "CC_10_PZ_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'CC_PZ', -100.)".format(stream, line),
		# "CC_10_PASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'CC_PASYM', -100.)".format(stream, line),
		# "CC_10_PTASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'CC_PTASYM', -100.)".format(stream, line),
		# "CC_10_PXASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'CC_PXASYM', -100.)".format(stream, line),
		# "CC_10_PYASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'CC_PYASYM', -100.)".format(stream, line),
		# "CC_10_PZASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'CC_PZASYM', -100.)".format(stream, line),
		# "CC_10_DELTAETA_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'CC_DELTAETA', -100.)".format(stream, line),
		# "CC_10_DELTAPHI_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'CC_DELTAPHI', -100.)".format(stream, line),
		# "CC_10_IT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'CC_IT', -100.)".format(stream, line),
		# "CC_10_MAXPT_Q_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'CC_MAXPT_Q', -100.)".format(stream, line),
		# "CC_10_MAXPT_PT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'CC_MAXPT_PT', -100.)".format(stream, line),
		# "CC_10_MAXPT_PX_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'CC_MAXPT_PX', -100.)".format(stream, line),
		# "CC_10_MAXPT_PY_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'CC_MAXPT_PY', -100.)".format(stream, line),
		# "CC_10_MAXPT_PZ_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'CC_MAXPT_PZ', -100.)".format(stream, line),
		# "CC_10_MAXPT_PE_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'CC_MAXPT_PE', -100.)".format(stream, line),
		# "NC_10_ANGLE_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'NC_ANGLE', -100.)".format(stream, line),
		# "NC_10_MULT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'NC_MULT', -100.)".format(stream, line),
		# "NC_10_SPT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'NC_SPT', -100.)".format(stream, line),
		# "NC_10_VPT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'NC_VPT', -100.)".format(stream, line),
		# "NC_10_PX_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'NC_PX', -100.)".format(stream, line),
		# "NC_10_PY_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'NC_PY', -100.)".format(stream, line),
		# "NC_10_PZ_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'NC_PZ', -100.)".format(stream, line),
		# "NC_10_PASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'NC_PASYM', -100.)".format(stream, line),
		# "NC_10_PTASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'NC_PTASYM', -100.)".format(stream, line),
		# "NC_10_PXASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'NC_PXASYM', -100.)".format(stream, line),
		# "NC_10_PYASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'NC_PYASYM', -100.)".format(stream, line),
		# "NC_10_PZASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'NC_PZASYM', -100.)".format(stream, line),
		# "NC_10_DELTAETA_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'NC_DELTAETA', -100.)".format(stream, line),
		# "NC_10_DELTAPHI_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'NC_DELTAPHI', -100.)".format(stream, line),
		# "NC_10_IT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'NC_IT', -100.)".format(stream, line),
		# "NC_10_MAXPT_PT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'NC_MAXPT_PT', -100.)".format(stream, line),
		# "NC_10_MAXPT_PX_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'NC_MAXPT_PX', -100.)".format(stream, line),
		# "NC_10_MAXPT_PY_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'NC_MAXPT_PY', -100.)".format(stream, line),
		# "NC_10_MAXPT_PZ_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'NC_MAXPT_PZ', -100.)".format(stream, line),
		# "CCNC_10_IT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone10', 'CCNC_IT', -100.)".format(stream, line),

		# # tau- (cone size = 1.0)
		# "CC_10_ANGLE_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'CC_ANGLE', -100.)".format(stream, line),
		# "CC_10_MULT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'CC_MULT', -100.)".format(stream, line),
		# "CC_10_SPT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'CC_SPT', -100.)".format(stream, line),
		# "CC_10_VPT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'CC_VPT', -100.)".format(stream, line),
		# "CC_10_PX_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'CC_PX', -100.)".format(stream, line),
		# "CC_10_PY_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'CC_PY', -100.)".format(stream, line),
		# "CC_10_PZ_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'CC_PZ', -100.)".format(stream, line),
		# "CC_10_PASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'CC_PASYM', -100.)".format(stream, line),
		# "CC_10_PTASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'CC_PTASYM', -100.)".format(stream, line),
		# "CC_10_PXASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'CC_PXASYM', -100.)".format(stream, line),
		# "CC_10_PYASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'CC_PYASYM', -100.)".format(stream, line),
		# "CC_10_PZASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'CC_PZASYM', -100.)".format(stream, line),
		# "CC_10_DELTAETA_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'CC_DELTAETA', -100.)".format(stream, line),
		# "CC_10_DELTAPHI_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'CC_DELTAPHI', -100.)".format(stream, line),
		# "CC_10_IT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'CC_IT', -100.)".format(stream, line),
		# "CC_10_MAXPT_Q_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'CC_MAXPT_Q', -100.)".format(stream, line),
		# "CC_10_MAXPT_PT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'CC_MAXPT_PT', -100.)".format(stream, line),
		# "CC_10_MAXPT_PX_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'CC_MAXPT_PX', -100.)".format(stream, line),
		# "CC_10_MAXPT_PY_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'CC_MAXPT_PY', -100.)".format(stream, line),
		# "CC_10_MAXPT_PZ_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'CC_MAXPT_PZ', -100.)".format(stream, line),
		# "CC_10_MAXPT_PE_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'CC_MAXPT_PE', -100.)".format(stream, line),
		# "NC_10_ANGLE_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'NC_ANGLE', -100.)".format(stream, line),
		# "NC_10_MULT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'NC_MULT', -100.)".format(stream, line),
		# "NC_10_SPT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'NC_SPT', -100.)".format(stream, line),
		# "NC_10_VPT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'NC_VPT', -100.)".format(stream, line),
		# "NC_10_PX_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'NC_PX', -100.)".format(stream, line),
		# "NC_10_PY_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'NC_PY', -100.)".format(stream, line),
		# "NC_10_PZ_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'NC_PZ', -100.)".format(stream, line),
		# "NC_10_PASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'NC_PASYM', -100.)".format(stream, line),
		# "NC_10_PTASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'NC_PTASYM', -100.)".format(stream, line),
		# "NC_10_PXASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'NC_PXASYM', -100.)".format(stream, line),
		# "NC_10_PYASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'NC_PYASYM', -100.)".format(stream, line),
		# "NC_10_PZASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'NC_PZASYM', -100.)".format(stream, line),
		# "NC_10_DELTAETA_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'NC_DELTAETA', -100.)".format(stream, line),
		# "NC_10_DELTAPHI_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'NC_DELTAPHI', -100.)".format(stream, line),
		# "NC_10_IT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'NC_IT', -100.)".format(stream, line),
		# "NC_10_MAXPT_PT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'NC_MAXPT_PT', -100.)".format(stream, line),
		# "NC_10_MAXPT_PX_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'NC_MAXPT_PX', -100.)".format(stream, line),
		# "NC_10_MAXPT_PY_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'NC_MAXPT_PY', -100.)".format(stream, line),
		# "NC_10_MAXPT_PZ_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'NC_MAXPT_PZ', -100.)".format(stream, line),
		# "CCNC_10_IT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone10', 'CCNC_IT', -100.)".format(stream, line),

		# # B+ (cone size = 1.5)
		# "CC_15_ANGLE_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_ANGLE', -100.)".format(stream, line),
		# "CC_15_MULT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_MULT', -100.)".format(stream, line),
		# "CC_15_SPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_SPT', -100.)".format(stream, line),
		# "CC_15_VPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_VPT', -100.)".format(stream, line),
		# "CC_15_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_PX', -100.)".format(stream, line),
		# "CC_15_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_PY', -100.)".format(stream, line),
		# "CC_15_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_PZ', -100.)".format(stream, line),
		# "CC_15_PASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_PASYM', -100.)".format(stream, line),
		# "CC_15_PTASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_PTASYM', -100.)".format(stream, line),
		# "CC_15_PXASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_PXASYM', -100.)".format(stream, line),
		# "CC_15_PYASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_PYASYM', -100.)".format(stream, line),
		# "CC_15_PZASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_PZASYM', -100.)".format(stream, line),
		# "CC_15_DELTAETA_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_DELTAETA', -100.)".format(stream, line),
		# "CC_15_DELTAPHI_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_DELTAPHI', -100.)".format(stream, line),
		# "CC_15_IT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_IT', -100.)".format(stream, line),
		# "CC_15_MAXPT_Q_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_MAXPT_Q', -100.)".format(stream, line),
		# "CC_15_MAXPT_PT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_MAXPT_PT', -100.)".format(stream, line),
		# "CC_15_MAXPT_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_MAXPT_PX', -100.)".format(stream, line),
		# "CC_15_MAXPT_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_MAXPT_PY', -100.)".format(stream, line),
		# "CC_15_MAXPT_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_MAXPT_PZ', -100.)".format(stream, line),
		# "CC_15_MAXPT_PE_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_MAXPT_PE', -100.)".format(stream, line),
		# "NC_15_ANGLE_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_ANGLE', -100.)".format(stream, line),
		# "NC_15_MULT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_MULT', -100.)".format(stream, line),
		# "NC_15_SPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_SPT', -100.)".format(stream, line),
		# "NC_15_VPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_VPT', -100.)".format(stream, line),
		# "NC_15_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_PX', -100.)".format(stream, line),
		# "NC_15_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_PY', -100.)".format(stream, line),
		# "NC_15_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_PZ', -100.)".format(stream, line),
		# "NC_15_PASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_PASYM', -100.)".format(stream, line),
		# "NC_15_PTASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_PTASYM', -100.)".format(stream, line),
		# "NC_15_PXASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_PXASYM', -100.)".format(stream, line),
		# "NC_15_PYASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_PYASYM', -100.)".format(stream, line),
		# "NC_15_PZASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_PZASYM', -100.)".format(stream, line),
		# "NC_15_DELTAETA_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_DELTAETA', -100.)".format(stream, line),
		# "NC_15_DELTAPHI_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_DELTAPHI', -100.)".format(stream, line),
		# "NC_15_IT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_IT', -100.)".format(stream, line),
		# "NC_15_MAXPT_PT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_MAXPT_PT', -100.)".format(stream, line),
		# "NC_15_MAXPT_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_MAXPT_PX', -100.)".format(stream, line),
		# "NC_15_MAXPT_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_MAXPT_PY', -100.)".format(stream, line),
		# "NC_15_MAXPT_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_MAXPT_PZ', -100.)".format(stream, line),
		# "CCNC_15_IT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CCNC_IT', -100.)".format(stream, line),

		# # K+ (cone size = 1.5)
		# "CC_15_ANGLE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_ANGLE', -100.)".format(stream, line),
		# "CC_15_MULT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_MULT', -100.)".format(stream, line),
		# "CC_15_SPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_SPT', -100.)".format(stream, line),
		# "CC_15_VPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_VPT', -100.)".format(stream, line),
		# "CC_15_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_PX', -100.)".format(stream, line),
		# "CC_15_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_PY', -100.)".format(stream, line),
		# "CC_15_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_PZ', -100.)".format(stream, line),
		# "CC_15_PASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_PASYM', -100.)".format(stream, line),
		# "CC_15_PTASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_PTASYM', -100.)".format(stream, line),
		# "CC_15_PXASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_PXASYM', -100.)".format(stream, line),
		# "CC_15_PYASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_PYASYM', -100.)".format(stream, line),
		# "CC_15_PZASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_PZASYM', -100.)".format(stream, line),
		# "CC_15_DELTAETA_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_DELTAETA', -100.)".format(stream, line),
		# "CC_15_DELTAPHI_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_DELTAPHI', -100.)".format(stream, line),
		# "CC_15_IT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_IT', -100.)".format(stream, line),
		# "CC_15_MAXPT_Q_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_MAXPT_Q', -100.)".format(stream, line),
		# "CC_15_MAXPT_PT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_MAXPT_PT', -100.)".format(stream, line),
		# "CC_15_MAXPT_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_MAXPT_PX', -100.)".format(stream, line),
		# "CC_15_MAXPT_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_MAXPT_PY', -100.)".format(stream, line),
		# "CC_15_MAXPT_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_MAXPT_PZ', -100.)".format(stream, line),
		# "CC_15_MAXPT_PE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_MAXPT_PE', -100.)".format(stream, line),
		# "NC_15_ANGLE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_ANGLE', -100.)".format(stream, line),
		# "NC_15_MULT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_MULT', -100.)".format(stream, line),
		# "NC_15_SPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_SPT', -100.)".format(stream, line),
		# "NC_15_VPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_VPT', -100.)".format(stream, line),
		# "NC_15_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_PX', -100.)".format(stream, line),
		# "NC_15_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_PY', -100.)".format(stream, line),
		# "NC_15_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_PZ', -100.)".format(stream, line),
		# "NC_15_PASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_PASYM', -100.)".format(stream, line),
		# "NC_15_PTASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_PTASYM', -100.)".format(stream, line),
		# "NC_15_PXASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_PXASYM', -100.)".format(stream, line),
		# "NC_15_PYASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_PYASYM', -100.)".format(stream, line),
		# "NC_15_PZASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_PZASYM', -100.)".format(stream, line),
		# "NC_15_DELTAETA_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_DELTAETA', -100.)".format(stream, line),
		# "NC_15_DELTAPHI_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_DELTAPHI', -100.)".format(stream, line),
		# "NC_15_IT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_IT', -100.)".format(stream, line),
		# "NC_15_MAXPT_PT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_MAXPT_PT', -100.)".format(stream, line),
		# "NC_15_MAXPT_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_MAXPT_PX', -100.)".format(stream, line),
		# "NC_15_MAXPT_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_MAXPT_PY', -100.)".format(stream, line),
		# "NC_15_MAXPT_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_MAXPT_PZ', -100.)".format(stream, line),
		# "CCNC_15_IT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CCNC_IT', -100.)".format(stream, line),

		# # tau+ (cone size = 1.5)
		# "CC_15_ANGLE_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'CC_ANGLE', -100.)".format(stream, line),
		# "CC_15_MULT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'CC_MULT', -100.)".format(stream, line),
		# "CC_15_SPT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'CC_SPT', -100.)".format(stream, line),
		# "CC_15_VPT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'CC_VPT', -100.)".format(stream, line),
		# "CC_15_PX_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'CC_PX', -100.)".format(stream, line),
		# "CC_15_PY_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'CC_PY', -100.)".format(stream, line),
		# "CC_15_PZ_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'CC_PZ', -100.)".format(stream, line),
		# "CC_15_PASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'CC_PASYM', -100.)".format(stream, line),
		# "CC_15_PTASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'CC_PTASYM', -100.)".format(stream, line),
		# "CC_15_PXASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'CC_PXASYM', -100.)".format(stream, line),
		# "CC_15_PYASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'CC_PYASYM', -100.)".format(stream, line),
		# "CC_15_PZASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'CC_PZASYM', -100.)".format(stream, line),
		# "CC_15_DELTAETA_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'CC_DELTAETA', -100.)".format(stream, line),
		# "CC_15_DELTAPHI_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'CC_DELTAPHI', -100.)".format(stream, line),
		# "CC_15_IT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'CC_IT', -100.)".format(stream, line),
		# "CC_15_MAXPT_Q_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'CC_MAXPT_Q', -100.)".format(stream, line),
		# "CC_15_MAXPT_PT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'CC_MAXPT_PT', -100.)".format(stream, line),
		# "CC_15_MAXPT_PX_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'CC_MAXPT_PX', -100.)".format(stream, line),
		# "CC_15_MAXPT_PY_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'CC_MAXPT_PY', -100.)".format(stream, line),
		# "CC_15_MAXPT_PZ_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'CC_MAXPT_PZ', -100.)".format(stream, line),
		# "CC_15_MAXPT_PE_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'CC_MAXPT_PE', -100.)".format(stream, line),
		# "NC_15_ANGLE_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'NC_ANGLE', -100.)".format(stream, line),
		# "NC_15_MULT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'NC_MULT', -100.)".format(stream, line),
		# "NC_15_SPT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'NC_SPT', -100.)".format(stream, line),
		# "NC_15_VPT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'NC_VPT', -100.)".format(stream, line),
		# "NC_15_PX_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'NC_PX', -100.)".format(stream, line),
		# "NC_15_PY_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'NC_PY', -100.)".format(stream, line),
		# "NC_15_PZ_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'NC_PZ', -100.)".format(stream, line),
		# "NC_15_PASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'NC_PASYM', -100.)".format(stream, line),
		# "NC_15_PTASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'NC_PTASYM', -100.)".format(stream, line),
		# "NC_15_PXASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'NC_PXASYM', -100.)".format(stream, line),
		# "NC_15_PYASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'NC_PYASYM', -100.)".format(stream, line),
		# "NC_15_PZASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'NC_PZASYM', -100.)".format(stream, line),
		# "NC_15_DELTAETA_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'NC_DELTAETA', -100.)".format(stream, line),
		# "NC_15_DELTAPHI_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'NC_DELTAPHI', -100.)".format(stream, line),
		# "NC_15_IT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'NC_IT', -100.)".format(stream, line),
		# "NC_15_MAXPT_PT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'NC_MAXPT_PT', -100.)".format(stream, line),
		# "NC_15_MAXPT_PX_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'NC_MAXPT_PX', -100.)".format(stream, line),
		# "NC_15_MAXPT_PY_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'NC_MAXPT_PY', -100.)".format(stream, line),
		# "NC_15_MAXPT_PZ_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'NC_MAXPT_PZ', -100.)".format(stream, line),
		# "CCNC_15_IT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone15', 'CCNC_IT', -100.)".format(stream, line),

		# # tau- (cone size = 1.5)
		# "CC_15_ANGLE_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'CC_ANGLE', -100.)".format(stream, line),
		# "CC_15_MULT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'CC_MULT', -100.)".format(stream, line),
		# "CC_15_SPT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'CC_SPT', -100.)".format(stream, line),
		# "CC_15_VPT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'CC_VPT', -100.)".format(stream, line),
		# "CC_15_PX_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'CC_PX', -100.)".format(stream, line),
		# "CC_15_PY_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'CC_PY', -100.)".format(stream, line),
		# "CC_15_PZ_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'CC_PZ', -100.)".format(stream, line),
		# "CC_15_PASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'CC_PASYM', -100.)".format(stream, line),
		# "CC_15_PTASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'CC_PTASYM', -100.)".format(stream, line),
		# "CC_15_PXASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'CC_PXASYM', -100.)".format(stream, line),
		# "CC_15_PYASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'CC_PYASYM', -100.)".format(stream, line),
		# "CC_15_PZASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'CC_PZASYM', -100.)".format(stream, line),
		# "CC_15_DELTAETA_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'CC_DELTAETA', -100.)".format(stream, line),
		# "CC_15_DELTAPHI_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'CC_DELTAPHI', -100.)".format(stream, line),
		# "CC_15_IT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'CC_IT', -100.)".format(stream, line),
		# "CC_15_MAXPT_Q_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'CC_MAXPT_Q', -100.)".format(stream, line),
		# "CC_15_MAXPT_PT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'CC_MAXPT_PT', -100.)".format(stream, line),
		# "CC_15_MAXPT_PX_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'CC_MAXPT_PX', -100.)".format(stream, line),
		# "CC_15_MAXPT_PY_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'CC_MAXPT_PY', -100.)".format(stream, line),
		# "CC_15_MAXPT_PZ_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'CC_MAXPT_PZ', -100.)".format(stream, line),
		# "CC_15_MAXPT_PE_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'CC_MAXPT_PE', -100.)".format(stream, line),
		# "NC_15_ANGLE_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'NC_ANGLE', -100.)".format(stream, line),
		# "NC_15_MULT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'NC_MULT', -100.)".format(stream, line),
		# "NC_15_SPT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'NC_SPT', -100.)".format(stream, line),
		# "NC_15_VPT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'NC_VPT', -100.)".format(stream, line),
		# "NC_15_PX_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'NC_PX', -100.)".format(stream, line),
		# "NC_15_PY_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'NC_PY', -100.)".format(stream, line),
		# "NC_15_PZ_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'NC_PZ', -100.)".format(stream, line),
		# "NC_15_PASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'NC_PASYM', -100.)".format(stream, line),
		# "NC_15_PTASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'NC_PTASYM', -100.)".format(stream, line),
		# "NC_15_PXASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'NC_PXASYM', -100.)".format(stream, line),
		# "NC_15_PYASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'NC_PYASYM', -100.)".format(stream, line),
		# "NC_15_PZASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'NC_PZASYM', -100.)".format(stream, line),
		# "NC_15_DELTAETA_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'NC_DELTAETA', -100.)".format(stream, line),
		# "NC_15_DELTAPHI_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'NC_DELTAPHI', -100.)".format(stream, line),
		# "NC_15_IT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'NC_IT', -100.)".format(stream, line),
		# "NC_15_MAXPT_PT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'NC_MAXPT_PT', -100.)".format(stream, line),
		# "NC_15_MAXPT_PX_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'NC_MAXPT_PX', -100.)".format(stream, line),
		# "NC_15_MAXPT_PY_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'NC_MAXPT_PY', -100.)".format(stream, line),
		# "NC_15_MAXPT_PZ_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'NC_MAXPT_PZ', -100.)".format(stream, line),
		# "CCNC_15_IT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone15', 'CCNC_IT', -100.)".format(stream, line),

		# # B+ (cone size = 2.0)
		# "CC_20_ANGLE_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_ANGLE', -100.)".format(stream, line),
		# "CC_20_MULT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_MULT', -100.)".format(stream, line),
		# "CC_20_SPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_SPT', -100.)".format(stream, line),
		# "CC_20_VPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_VPT', -100.)".format(stream, line),
		# "CC_20_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_PX', -100.)".format(stream, line),
		# "CC_20_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_PY', -100.)".format(stream, line),
		# "CC_20_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_PZ', -100.)".format(stream, line),
		# "CC_20_PASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_PASYM', -100.)".format(stream, line),
		# "CC_20_PTASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_PTASYM', -100.)".format(stream, line),
		# "CC_20_PXASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_PXASYM', -100.)".format(stream, line),
		# "CC_20_PYASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_PYASYM', -100.)".format(stream, line),
		# "CC_20_PZASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_PZASYM', -100.)".format(stream, line),
		# "CC_20_DELTAETA_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_DELTAETA', -100.)".format(stream, line),
		# "CC_20_DELTAPHI_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_DELTAPHI', -100.)".format(stream, line),
		# "CC_20_IT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_IT', -100.)".format(stream, line),
		# "CC_20_MAXPT_Q_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_MAXPT_Q', -100.)".format(stream, line),
		# "CC_20_MAXPT_PT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_MAXPT_PT', -100.)".format(stream, line),
		# "CC_20_MAXPT_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_MAXPT_PX', -100.)".format(stream, line),
		# "CC_20_MAXPT_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_MAXPT_PY', -100.)".format(stream, line),
		# "CC_20_MAXPT_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_MAXPT_PZ', -100.)".format(stream, line),
		# "CC_20_MAXPT_PE_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_MAXPT_PE', -100.)".format(stream, line),
		# "NC_20_ANGLE_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_ANGLE', -100.)".format(stream, line),
		# "NC_20_MULT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_MULT', -100.)".format(stream, line),
		# "NC_20_SPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_SPT', -100.)".format(stream, line),
		# "NC_20_VPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_VPT', -100.)".format(stream, line),
		# "NC_20_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_PX', -100.)".format(stream, line),
		# "NC_20_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_PY', -100.)".format(stream, line),
		# "NC_20_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_PZ', -100.)".format(stream, line),
		# "NC_20_PASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_PASYM', -100.)".format(stream, line),
		# "NC_20_PTASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_PTASYM', -100.)".format(stream, line),
		# "NC_20_PXASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_PXASYM', -100.)".format(stream, line),
		# "NC_20_PYASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_PYASYM', -100.)".format(stream, line),
		# "NC_20_PZASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_PZASYM', -100.)".format(stream, line),
		# "NC_20_DELTAETA_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_DELTAETA', -100.)".format(stream, line),
		# "NC_20_DELTAPHI_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_DELTAPHI', -100.)".format(stream, line),
		# "NC_20_IT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_IT', -100.)".format(stream, line),
		# "NC_20_MAXPT_PT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_MAXPT_PT', -100.)".format(stream, line),
		# "NC_20_MAXPT_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_MAXPT_PX', -100.)".format(stream, line),
		# "NC_20_MAXPT_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_MAXPT_PY', -100.)".format(stream, line),
		# "NC_20_MAXPT_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_MAXPT_PZ', -100.)".format(stream, line),
		# "CCNC_20_IT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CCNC_IT', -100.)".format(stream, line),

		# # K+ (cone size = 2.0)
		# "CC_20_ANGLE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_ANGLE', -100.)".format(stream, line),
		# "CC_20_MULT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_MULT', -100.)".format(stream, line),
		# "CC_20_SPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_SPT', -100.)".format(stream, line),
		# "CC_20_VPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_VPT', -100.)".format(stream, line),
		# "CC_20_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_PX', -100.)".format(stream, line),
		# "CC_20_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_PY', -100.)".format(stream, line),
		# "CC_20_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_PZ', -100.)".format(stream, line),
		# "CC_20_PASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_PASYM', -100.)".format(stream, line),
		# "CC_20_PTASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_PTASYM', -100.)".format(stream, line),
		# "CC_20_PXASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_PXASYM', -100.)".format(stream, line),
		# "CC_20_PYASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_PYASYM', -100.)".format(stream, line),
		# "CC_20_PZASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_PZASYM', -100.)".format(stream, line),
		# "CC_20_DELTAETA_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_DELTAETA', -100.)".format(stream, line),
		# "CC_20_DELTAPHI_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_DELTAPHI', -100.)".format(stream, line),
		# "CC_20_IT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_IT', -100.)".format(stream, line),
		# "CC_20_MAXPT_Q_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_MAXPT_Q', -100.)".format(stream, line),
		# "CC_20_MAXPT_PT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_MAXPT_PT', -100.)".format(stream, line),
		# "CC_20_MAXPT_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_MAXPT_PX', -100.)".format(stream, line),
		# "CC_20_MAXPT_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_MAXPT_PY', -100.)".format(stream, line),
		# "CC_20_MAXPT_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_MAXPT_PZ', -100.)".format(stream, line),
		# "CC_20_MAXPT_PE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_MAXPT_PE', -100.)".format(stream, line),
		# "NC_20_ANGLE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_ANGLE', -100.)".format(stream, line),
		# "NC_20_MULT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_MULT', -100.)".format(stream, line),
		# "NC_20_SPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_SPT', -100.)".format(stream, line),
		# "NC_20_VPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_VPT', -100.)".format(stream, line),
		# "NC_20_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_PX', -100.)".format(stream, line),
		# "NC_20_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_PY', -100.)".format(stream, line),
		# "NC_20_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_PZ', -100.)".format(stream, line),
		# "NC_20_PASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_PASYM', -100.)".format(stream, line),
		# "NC_20_PTASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_PTASYM', -100.)".format(stream, line),
		# "NC_20_PXASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_PXASYM', -100.)".format(stream, line),
		# "NC_20_PYASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_PYASYM', -100.)".format(stream, line),
		# "NC_20_PZASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_PZASYM', -100.)".format(stream, line),
		# "NC_20_DELTAETA_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_DELTAETA', -100.)".format(stream, line),
		# "NC_20_DELTAPHI_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_DELTAPHI', -100.)".format(stream, line),
		# "NC_20_IT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_IT', -100.)".format(stream, line),
		# "NC_20_MAXPT_PT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_MAXPT_PT', -100.)".format(stream, line),
		# "NC_20_MAXPT_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_MAXPT_PX', -100.)".format(stream, line),
		# "NC_20_MAXPT_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_MAXPT_PY', -100.)".format(stream, line),
		# "NC_20_MAXPT_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_MAXPT_PZ', -100.)".format(stream, line),
		# "CCNC_20_IT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CCNC_IT', -100.)".format(stream, line),

		# # tau+ (cone size = 2.0)
		# "CC_20_ANGLE_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'CC_ANGLE', -100.)".format(stream, line),
		# "CC_20_MULT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'CC_MULT', -100.)".format(stream, line),
		# "CC_20_SPT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'CC_SPT', -100.)".format(stream, line),
		# "CC_20_VPT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'CC_VPT', -100.)".format(stream, line),
		# "CC_20_PX_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'CC_PX', -100.)".format(stream, line),
		# "CC_20_PY_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'CC_PY', -100.)".format(stream, line),
		# "CC_20_PZ_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'CC_PZ', -100.)".format(stream, line),
		# "CC_20_PASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'CC_PASYM', -100.)".format(stream, line),
		# "CC_20_PTASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'CC_PTASYM', -100.)".format(stream, line),
		# "CC_20_PXASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'CC_PXASYM', -100.)".format(stream, line),
		# "CC_20_PYASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'CC_PYASYM', -100.)".format(stream, line),
		# "CC_20_PZASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'CC_PZASYM', -100.)".format(stream, line),
		# "CC_20_DELTAETA_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'CC_DELTAETA', -100.)".format(stream, line),
		# "CC_20_DELTAPHI_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'CC_DELTAPHI', -100.)".format(stream, line),
		# "CC_20_IT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'CC_IT', -100.)".format(stream, line),
		# "CC_20_MAXPT_Q_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'CC_MAXPT_Q', -100.)".format(stream, line),
		# "CC_20_MAXPT_PT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'CC_MAXPT_PT', -100.)".format(stream, line),
		# "CC_20_MAXPT_PX_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'CC_MAXPT_PX', -100.)".format(stream, line),
		# "CC_20_MAXPT_PY_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'CC_MAXPT_PY', -100.)".format(stream, line),
		# "CC_20_MAXPT_PZ_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'CC_MAXPT_PZ', -100.)".format(stream, line),
		# "CC_20_MAXPT_PE_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'CC_MAXPT_PE', -100.)".format(stream, line),
		# "NC_20_ANGLE_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'NC_ANGLE', -100.)".format(stream, line),
		# "NC_20_MULT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'NC_MULT', -100.)".format(stream, line),
		# "NC_20_SPT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'NC_SPT', -100.)".format(stream, line),
		# "NC_20_VPT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'NC_VPT', -100.)".format(stream, line),
		# "NC_20_PX_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'NC_PX', -100.)".format(stream, line),
		# "NC_20_PY_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'NC_PY', -100.)".format(stream, line),
		# "NC_20_PZ_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'NC_PZ', -100.)".format(stream, line),
		# "NC_20_PASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'NC_PASYM', -100.)".format(stream, line),
		# "NC_20_PTASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'NC_PTASYM', -100.)".format(stream, line),
		# "NC_20_PXASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'NC_PXASYM', -100.)".format(stream, line),
		# "NC_20_PYASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'NC_PYASYM', -100.)".format(stream, line),
		# "NC_20_PZASYM_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'NC_PZASYM', -100.)".format(stream, line),
		# "NC_20_DELTAETA_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'NC_DELTAETA', -100.)".format(stream, line),
		# "NC_20_DELTAPHI_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'NC_DELTAPHI', -100.)".format(stream, line),
		# "NC_20_IT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'NC_IT', -100.)".format(stream, line),
		# "NC_20_MAXPT_PT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'NC_MAXPT_PT', -100.)".format(stream, line),
		# "NC_20_MAXPT_PX_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'NC_MAXPT_PX', -100.)".format(stream, line),
		# "NC_20_MAXPT_PY_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'NC_MAXPT_PY', -100.)".format(stream, line),
		# "NC_20_MAXPT_PZ_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'NC_MAXPT_PZ', -100.)".format(stream, line),
		# "CCNC_20_IT_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_ConeIsoInfo_Cone20', 'CCNC_IT', -100.)".format(stream, line),

		# # tau- (cone size = 2.0)
		# "CC_20_ANGLE_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'CC_ANGLE', -100.)".format(stream, line),
		# "CC_20_MULT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'CC_MULT', -100.)".format(stream, line),
		# "CC_20_SPT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'CC_SPT', -100.)".format(stream, line),
		# "CC_20_VPT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'CC_VPT', -100.)".format(stream, line),
		# "CC_20_PX_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'CC_PX', -100.)".format(stream, line),
		# "CC_20_PY_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'CC_PY', -100.)".format(stream, line),
		# "CC_20_PZ_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'CC_PZ', -100.)".format(stream, line),
		# "CC_20_PASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'CC_PASYM', -100.)".format(stream, line),
		# "CC_20_PTASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'CC_PTASYM', -100.)".format(stream, line),
		# "CC_20_PXASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'CC_PXASYM', -100.)".format(stream, line),
		# "CC_20_PYASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'CC_PYASYM', -100.)".format(stream, line),
		# "CC_20_PZASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'CC_PZASYM', -100.)".format(stream, line),
		# "CC_20_DELTAETA_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'CC_DELTAETA', -100.)".format(stream, line),
		# "CC_20_DELTAPHI_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'CC_DELTAPHI', -100.)".format(stream, line),
		# "CC_20_IT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'CC_IT', -100.)".format(stream, line),
		# "CC_20_MAXPT_Q_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'CC_MAXPT_Q', -100.)".format(stream, line),
		# "CC_20_MAXPT_PT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'CC_MAXPT_PT', -100.)".format(stream, line),
		# "CC_20_MAXPT_PX_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'CC_MAXPT_PX', -100.)".format(stream, line),
		# "CC_20_MAXPT_PY_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'CC_MAXPT_PY', -100.)".format(stream, line),
		# "CC_20_MAXPT_PZ_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'CC_MAXPT_PZ', -100.)".format(stream, line),
		# "CC_20_MAXPT_PE_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'CC_MAXPT_PE', -100.)".format(stream, line),
		# "NC_20_ANGLE_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'NC_ANGLE', -100.)".format(stream, line),
		# "NC_20_MULT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'NC_MULT', -100.)".format(stream, line),
		# "NC_20_SPT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'NC_SPT', -100.)".format(stream, line),
		# "NC_20_VPT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'NC_VPT', -100.)".format(stream, line),
		# "NC_20_PX_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'NC_PX', -100.)".format(stream, line),
		# "NC_20_PY_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'NC_PY', -100.)".format(stream, line),
		# "NC_20_PZ_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'NC_PZ', -100.)".format(stream, line),
		# "NC_20_PASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'NC_PASYM', -100.)".format(stream, line),
		# "NC_20_PTASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'NC_PTASYM', -100.)".format(stream, line),
		# "NC_20_PXASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'NC_PXASYM', -100.)".format(stream, line),
		# "NC_20_PYASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'NC_PYASYM', -100.)".format(stream, line),
		# "NC_20_PZASYM_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'NC_PZASYM', -100.)".format(stream, line),
		# "NC_20_DELTAETA_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'NC_DELTAETA', -100.)".format(stream, line),
		# "NC_20_DELTAPHI_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'NC_DELTAPHI', -100.)".format(stream, line),
		# "NC_20_IT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'NC_IT', -100.)".format(stream, line),
		# "NC_20_MAXPT_PT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'NC_MAXPT_PT', -100.)".format(stream, line),
		# "NC_20_MAXPT_PX_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'NC_MAXPT_PX', -100.)".format(stream, line),
		# "NC_20_MAXPT_PY_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'NC_MAXPT_PY', -100.)".format(stream, line),
		# "NC_20_MAXPT_PZ_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'NC_MAXPT_PZ', -100.)".format(stream, line),
		# "CCNC_20_IT_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_ConeIsoInfo_Cone20', 'CCNC_IT', -100.)".format(stream, line),

		# Vertex isolation BDT
		# B+
		"VTXISOBDTHARDFIRSTVALUE_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_VertexIsoBDTInfo', 'VTXISOBDTHARDFIRSTVALUE', -100.)".format(stream, line),
		"VTXISOBDTHARDSECONDVALUE_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_VertexIsoBDTInfo', 'VTXISOBDTHARDSECONDVALUE', -100.)".format(stream, line),
		"VTXISOBDTHARDTHIRDVALUE_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_VertexIsoBDTInfo', 'VTXISOBDTHARDTHIRDVALUE', -100.)".format(stream, line),

		# K+
		"VTXISOBDTHARDFIRSTVALUE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_VertexIsoBDTInfo', 'VTXISOBDTHARDFIRSTVALUE', -100.)".format(stream, line),
		"VTXISOBDTHARDSECONDVALUE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_VertexIsoBDTInfo', 'VTXISOBDTHARDSECONDVALUE', -100.)".format(stream, line),
		"VTXISOBDTHARDTHIRDVALUE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_VertexIsoBDTInfo', 'VTXISOBDTHARDTHIRDVALUE', -100.)".format(stream, line),

		# tau+
		"VTXISOBDTHARDFIRSTVALUE_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_VertexIsoBDTInfo', 'VTXISOBDTHARDFIRSTVALUE', -100.)".format(stream, line),
		"VTXISOBDTHARDSECONDVALUE_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_VertexIsoBDTInfo', 'VTXISOBDTHARDSECONDVALUE', -100.)".format(stream, line),
		"VTXISOBDTHARDTHIRDVALUE_taup" : "RELINFO('/Event/{0}/Phys/{1}/Taup_VertexIsoBDTInfo', 'VTXISOBDTHARDTHIRDVALUE', -100.)".format(stream, line),

		# tau-
		"VTXISOBDTHARDFIRSTVALUE_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_VertexIsoBDTInfo', 'VTXISOBDTHARDFIRSTVALUE', -100.)".format(stream, line),
		"VTXISOBDTHARDSECONDVALUE_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_VertexIsoBDTInfo', 'VTXISOBDTHARDSECONDVALUE', -100.)".format(stream, line),
		"VTXISOBDTHARDTHIRDVALUE_taum" : "RELINFO('/Event/{0}/Phys/{1}/Taum_VertexIsoBDTInfo', 'VTXISOBDTHARDTHIRDVALUE', -100.)".format(stream, line),

		# Track isolation BDT
		# K+
		"TRKISOBDTFIRSTVALUE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_TrackIsoBDTInfo', 'TRKISOBDTFIRSTVALUE', -100.)".format(stream, line),
		"TRKISOBDTSECONDVALUE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_TrackIsoBDTInfo', 'TRKISOBDTSECONDVALUE', -100.)".format(stream, line),
		"TRKISOBDTTHIRDVALUE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_TrackIsoBDTInfo', 'TRKISOBDTTHIRDVALUE', -100.)".format(stream, line),

		# taup pi1
		"TRKISOBDTFIRSTVALUE_taup_pi1" : "RELINFO('/Event/{0}/Phys/{1}/Taup_pi1_TrackIsoBDTInfo', 'TRKISOBDTFIRSTVALUE', -100.)".format(stream, line),
		"TRKISOBDTSECONDVALUE_taup_pi1" : "RELINFO('/Event/{0}/Phys/{1}/Taup_pi1_TrackIsoBDTInfo', 'TRKISOBDTSECONDVALUE', -100.)".format(stream, line),
		"TRKISOBDTTHIRDVALUE_taup_pi1" : "RELINFO('/Event/{0}/Phys/{1}/Taup_pi1_TrackIsoBDTInfo', 'TRKISOBDTTHIRDVALUE', -100.)".format(stream, line),

		# taup pi2
		"TRKISOBDTFIRSTVALUE_taup_pi2" : "RELINFO('/Event/{0}/Phys/{1}/Taup_pi2_TrackIsoBDTInfo', 'TRKISOBDTFIRSTVALUE', -100.)".format(stream, line),
		"TRKISOBDTSECONDVALUE_taup_pi2" : "RELINFO('/Event/{0}/Phys/{1}/Taup_pi2_TrackIsoBDTInfo', 'TRKISOBDTSECONDVALUE', -100.)".format(stream, line),
		"TRKISOBDTTHIRDVALUE_taup_pi2" : "RELINFO('/Event/{0}/Phys/{1}/Taup_pi2_TrackIsoBDTInfo', 'TRKISOBDTTHIRDVALUE', -100.)".format(stream, line),

		# taup pi3
		"TRKISOBDTFIRSTVALUE_taup_pi3" : "RELINFO('/Event/{0}/Phys/{1}/Taup_pi3_TrackIsoBDTInfo', 'TRKISOBDTFIRSTVALUE', -100.)".format(stream, line),
		"TRKISOBDTSECONDVALUE_taup_pi3" : "RELINFO('/Event/{0}/Phys/{1}/Taup_pi3_TrackIsoBDTInfo', 'TRKISOBDTSECONDVALUE', -100.)".format(stream, line),
		"TRKISOBDTTHIRDVALUE_taup_pi3" : "RELINFO('/Event/{0}/Phys/{1}/Taup_pi3_TrackIsoBDTInfo', 'TRKISOBDTTHIRDVALUE', -100.)".format(stream, line),

		# taum pi1
		"TRKISOBDTFIRSTVALUE_taum_pi1" : "RELINFO('/Event/{0}/Phys/{1}/Taum_pi1_TrackIsoBDTInfo', 'TRKISOBDTFIRSTVALUE', -100.)".format(stream, line),
		"TRKISOBDTSECONDVALUE_taum_pi1" : "RELINFO('/Event/{0}/Phys/{1}/Taum_pi1_TrackIsoBDTInfo', 'TRKISOBDTSECONDVALUE', -100.)".format(stream, line),
		"TRKISOBDTTHIRDVALUE_taum_pi1" : "RELINFO('/Event/{0}/Phys/{1}/Taum_pi1_TrackIsoBDTInfo', 'TRKISOBDTTHIRDVALUE', -100.)".format(stream, line),

		# taum pi2
		"TRKISOBDTFIRSTVALUE_taum_pi2" : "RELINFO('/Event/{0}/Phys/{1}/Taum_pi2_TrackIsoBDTInfo', 'TRKISOBDTFIRSTVALUE', -100.)".format(stream, line),
		"TRKISOBDTSECONDVALUE_taum_pi2" : "RELINFO('/Event/{0}/Phys/{1}/Taum_pi2_TrackIsoBDTInfo', 'TRKISOBDTSECONDVALUE', -100.)".format(stream, line),
		"TRKISOBDTTHIRDVALUE_taum_pi2" : "RELINFO('/Event/{0}/Phys/{1}/Taum_pi2_TrackIsoBDTInfo', 'TRKISOBDTTHIRDVALUE', -100.)".format(stream, line),

		# taum pi3
		"TRKISOBDTFIRSTVALUE_taum_pi3" : "RELINFO('/Event/{0}/Phys/{1}/Taum_pi3_TrackIsoBDTInfo', 'TRKISOBDTFIRSTVALUE', -100.)".format(stream, line),
		"TRKISOBDTSECONDVALUE_taum_pi3" : "RELINFO('/Event/{0}/Phys/{1}/Taum_pi3_TrackIsoBDTInfo', 'TRKISOBDTSECONDVALUE', -100.)".format(stream, line),
		"TRKISOBDTTHIRDVALUE_taum_pi3" : "RELINFO('/Event/{0}/Phys/{1}/Taum_pi3_TrackIsoBDTInfo', 'TRKISOBDTTHIRDVALUE', -100.)".format(stream, line),

		# B2Ksttautau tau isolation BDT
		"B2Ksttautau_ISOBDTFIRSTVALUE_taup" : "RELINFO('/Event/{0}/Phys/{1}/B2KstTauTau_TauIsolationBDT', 'BKSTTAUTAUTAUISOBDTFIRSTVALUETAUP', -100.)".format(stream, line),
		"B2Ksttautau_ISOBDTSECONDVALUE_taup" : "RELINFO('/Event/{0}/Phys/{1}/B2KstTauTau_TauIsolationBDT', 'BKSTTAUTAUTAUISOBDTSECONDVALUETAUP', -100.)".format(stream, line),
		"B2Ksttautau_ISOBDTTHIRDVALUE_taup" : "RELINFO('/Event/{0}/Phys/{1}/B2KstTauTau_TauIsolationBDT', 'BKSTTAUTAUTAUISOBDTTHIRDVALUETAUP', -100.)".format(stream, line),
		"B2Ksttautau_ISOBDTFIRSTVALUE_taum" : "RELINFO('/Event/{0}/Phys/{1}/B2KstTauTau_TauIsolationBDT', 'BKSTTAUTAUTAUISOBDTFIRSTVALUETAUM', -100.)".format(stream, line),
		"B2Ksttautau_ISOBDTSECONDVALUE_taum" : "RELINFO('/Event/{0}/Phys/{1}/B2KstTauTau_TauIsolationBDT', 'BKSTTAUTAUTAUISOBDTSECONDVALUETAUM', -100.)".format(stream, line),
		"B2Ksttautau_ISOBDTTHIRDVALUE_taum" : "RELINFO('/Event/{0}/Phys/{1}/B2KstTauTau_TauIsolationBDT', 'BKSTTAUTAUTAUISOBDTTHIRDVALUETAUM', -100.)".format(stream, line)
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
    	"VCHI2PDOF":"VFASPF(VCHI2PDOF)"
        # "AM" : "PFUNA(AM)",
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

# DecayTreeFitter / standalone fitter input branches
def addDTF(branch):
	# my decay fit (inputs for DECAY_FIT.C)
	df = branch.addTupleTool("TupleToolDecayFit/df")
	df.constrainToOriginVertex = True

	# # modified DTF (MLP init)
	# dtf_mlp = branch.addTupleTool("TupleToolKTauTauDTF/dtf_mlp")
	# dtf_mlp.daughtersToConstrain = ['nu_tau', 'tau+']
	# dtf_mlp.UseFullTreeInName = True
	# dtf_mlp.UpdateDaughters = True
	# dtf_mlp.Verbose = True
	# dtf_mlp.constrainToOriginVertex = True
	# dtf_mlp.init = 0

	# modified DTF (Low chi2 init)
	# dtf_low = branch.addTupleTool("TupleToolKTauTauDTF/dtf")
	# dtf_low.daughtersToConstrain = ['nu_tau', 'tau+']
	# dtf_low.UseFullTreeInName = True
	# dtf_low.UpdateDaughters = True
	# dtf_low.Verbose = True
	# dtf_low.constrainToOriginVertex = True
	# dtf_low.init = 1

def addTISTOS( branch ):
    ttb = branch.addTupleTool("TupleToolTISTOS")
    ttb.VerboseL0=True # if you want to fill info for each trigger line
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
    addDTF(DecayTreeTuple.Bp) 
    #add some helpful variables for B0
    addTopInfo(DecayTreeTuple.Bp)
    #add isolation info
    addIsoInfo(DecayTreeTuple, stream, line, year)
    #add submasses
    #addSubmassInfo(DecayTreeTuple.Bp)
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
    if isMC:
        addTupleToolL0Calo(DecayTreeTuple.taup_pi1)
        addTupleToolL0Calo(DecayTreeTuple.taup_pi2)
        addTupleToolL0Calo(DecayTreeTuple.taup_pi3)
        addTupleToolL0Calo(DecayTreeTuple.taum_pi1)
        addTupleToolL0Calo(DecayTreeTuple.taum_pi2)
        addTupleToolL0Calo(DecayTreeTuple.taum_pi3)
        addTupleToolL0Calo(DecayTreeTuple.Kp)

# Configure DaVinci
if isMC:
	# DaVinci().TupleFile = 'DVntuple_MC_BuDDKp_cocktail_'+year+'_'+pol+'.root'
	DaVinci().TupleFile = '/panfs/felician/BuDDKp_cocktail/DVntuple_MC_'+year+'_'+pol+'.root'
else:
	DaVinci().TupleFile = 'DVntuple_data_BuDDKp_cocktail_'+year+'_'+pol+'.root'
	# DaVinci().TupleFile = '/panfs/felician/BuDDKp_cocktail//DVntuple_data_'+year+'_'+pol+'.root'

# Stream and stripping line we want to use
if isMC:
	stream = 'B2KTauTau.Strip' # for filtered MC
else:	
	stream = 'Bhadron'
line = 'B2KTauTauLine'
line_SS = 'B2KTauTauSSLine'

if isMC:
	DaVinci().Simulation = True
	DaVinci().InputType = 'DST'

	# CondDBtag / DDDBtag specify the exact detector conditions with which the MC was generated. They are specified in the downloaded DST file.
	if year == '2016':
		DaVinci().DDDBtag = "dddb-20220927-2016"
		if pol == 'MagUp':
			DaVinci().CondDBtag = "sim-20201113-6-vc-mu100-Sim10" 
		elif pol == 'MagDown':
			DaVinci().CondDBtag = "sim-20201113-6-vc-md100-Sim10"
	elif year == '2017':
		DaVinci().DDDBtag = "dddb-20220927-2017"
		if pol == 'MagUp':
			DaVinci().CondDBtag = "sim-20201113-7-vc-mu100-Sim10"
		elif pol == 'MagDown':
			DaVinci().CondDBtag = "sim-20201113-7-vc-md100-Sim10"
	elif year == '2018':
		DaVinci().DDDBtag = "dddb-20220927-2018"
		if pol == 'MagUp':
			DaVinci().CondDBtag = "sim-20201113-8-vc-mu100-Sim10" 
		elif pol == 'MagDown':
			DaVinci().CondDBtag = "sim-20201113-8-vc-md100-Sim10"
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
dtt_SS = DecayTreeTuple('ntuple_SS')
mcdtt = MCDecayTreeTuple('mc_ntuple')
# B+ -> D0D0 K+
mcdtt_BuD0D0Kp = MCDecayTreeTuple('mc_ntuple_BuD0D0Kp')
mcdtt_BuD0starD0Kp = MCDecayTreeTuple('mc_ntuple_BuD0starD0Kp')
mcdtt_BuD0D0starKp = MCDecayTreeTuple('mc_ntuple_BuD0D0starKp')
mcdtt_BuD0starD0starKp = MCDecayTreeTuple('mc_ntuple_BuD0starD0starKp')
# B+ -> D+D- K+
mcdtt_BuDpDmKp = MCDecayTreeTuple('mc_ntuple_BuDpDmKp')
mcdtt_BuDpstarDmKp = MCDecayTreeTuple('mc_ntuple_BuDpstarDmKp')
mcdtt_BuDpDmstarKp = MCDecayTreeTuple('mc_ntuple_BuDpDmstarKp')
mcdtt_BuDpstarDmstarKp = MCDecayTreeTuple('mc_ntuple_BuDpstarDmstarKp')
# B+ -> Ds+Ds- K+
mcdtt_BuDsDsKp = MCDecayTreeTuple('mc_ntuple_BuDsDsKp')
mcdtt_BuDsstarDsKp = MCDecayTreeTuple('mc_ntuple_BuDsstarDsKp')
mcdtt_BuDsDsstarKp = MCDecayTreeTuple('mc_ntuple_BuDsDsstarKp')
mcdtt_BuDsstarDsstarKp = MCDecayTreeTuple('mc_ntuple_BuDsstarDsstarKp')

# Preselections to reduce the size of the tuples (in data)
filtercode = "(NINGENERATION( (M < 750) & (ABSID=='tau+') , 1) == 0) & (NINGENERATION( (M > 1650) & (ABSID=='tau+') , 1) == 0) & (NINGENERATION( (BPVVD < 4) & (ABSID=='B+') , 0) == 0) & ( ((RELINFO('/Event/{0}/Phys/{1}/BVars_VertexIsoInfo', 'VTXISODCHI2MASSONETRACK', -100) > 3600) & (RELINFO('/Event/{0}/Phys/{1}/BVars_VertexIsoBDTInfo', 'VTXISOBDTHARDFIRSTVALUE', -100.) < 0) ) | ((RELINFO('/Event/{0}/Phys/{2}/BVars_VertexIsoInfo', 'VTXISODCHI2MASSONETRACK', -100) > 3600) & (RELINFO('/Event/{0}/Phys/{2}/BVars_VertexIsoBDTInfo', 'VTXISOBDTHARDFIRSTVALUE', -100.) < 0) ) )".format(stream,line,line_SS)

if applySel:
	# if isMC:
	# 	b_strip = AutomaticData('/Event/{0}/Phys/{1}/Particles'.format(stream, line))
	# if not isMC:
	b_strip    = AutomaticData('/Phys/{0}/Particles'.format(line)) 
	b_strip_SS = AutomaticData('/Phys/{0}/Particles'.format(line_SS)) 

	b_Filter    = FilterSelection("b_Filter", b_strip, Code = filtercode)
	dtt.Inputs = [b_Filter.outputLocation()]

	# if not isMC:
	b_Filter_SS = FilterSelection("b_Filter_SS", b_strip_SS, Code = filtercode)
	dtt_SS.Inputs = [b_Filter_SS.outputLocation()]

# Do not apply rectangular cuts to MC @ DaVinci stage:
if not applySel:
	if isMC:
		dtt.Inputs = ['/Event/{0}/Phys/{1}/Particles'.format(stream, line)] 
	if not isMC:
		dtt.Inputs = ['/Phys/{0}/Particles'.format(line)]
		dtt_SS.Inputs = ['/Phys/{0}/Particles'.format(line_SS)]

dtt.setDescriptorTemplate('${Bp}[B+ -> ${taup}(tau+ -> ${taup_pi1}pi+ ${taup_pi2}pi- ${taup_pi3}pi+) ${taum}(tau- -> ${taum_pi1}pi- ${taum_pi2}pi+ ${taum_pi3}pi-) ${Kp}K+]CC')
dtt_SS.setDescriptorTemplate('(${Bp}[B+ -> ${taup}(tau+ -> ${taup_pi1}pi+ ${taup_pi2}pi- ${taup_pi3}pi+) ${taum}(tau+ -> ${taum_pi1}pi+ ${taum_pi2}pi- ${taum_pi3}pi+) ${Kp}K+]CC)'
						  '|| (${Bp}[B+ -> ${taup}(tau+ -> ${taup_pi1}pi+ ${taup_pi2}pi- ${taup_pi3}pi+) ${taum}(tau+ -> ${taum_pi1}pi+ ${taum_pi2}pi- ${taum_pi3}pi+) ${Kp}K-]CC)')
# Generated 
mcdtt.setDescriptorTemplate('${Bp}[B+ => ${D1}(Charm ==> ${D1_pi1}pi+ ${D1_pi2}pi- ${D1_pi3}pi+ {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${D2}(Charm ==> ${D2_pi1}pi- ${D2_pi2}pi+ ${D2_pi3}pi- {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${Kp}K+]CC')

# B+ -> D0D0 K+
mcdtt_BuD0D0Kp.setDescriptorTemplate('(${Bp}[B+ => ${D1}(D~0 ==> ${D1_pi1}pi+ ${D1_pi2}pi- ${D1_pi3}pi+ {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${D2}(D0 ==> ${D2_pi1}pi- ${D2_pi2}pi+ ${D2_pi3}pi- {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${Kp}K+]CC) || (${Bp}[B+ => ${D1}(D~0 ==> ${D1_pi1}pi- ${D1_pi2}pi+ ${D1_pi3}pi- {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${D2}(D0 ==> ${D2_pi1}pi+ ${D2_pi2}pi- ${D2_pi3}pi+ {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${Kp}K+]CC)')
mcdtt_BuD0starD0Kp.setDescriptorTemplate('(${Bp}[B+ => ${D1}(D*(2007)~0 ==> ${D1_pi1}pi+ ${D1_pi2}pi- ${D1_pi3}pi+ {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${D2}(D0 ==> ${D2_pi1}pi- ${D2_pi2}pi+ ${D2_pi3}pi- {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${Kp}K+]CC) || (${Bp}[B+ => ${D1}(D*(2007)~0 ==> ${D1_pi1}pi- ${D1_pi2}pi+ ${D1_pi3}pi- {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${D2}(D0 ==> ${D2_pi1}pi+ ${D2_pi2}pi- ${D2_pi3}pi+ {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${Kp}K+]CC)')
mcdtt_BuD0D0starKp.setDescriptorTemplate('(${Bp}[B+ => ${D1}(D~0 ==> ${D1_pi1}pi+ ${D1_pi2}pi- ${D1_pi3}pi+ {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${D2}(D*(2007)0 ==> ${D2_pi1}pi- ${D2_pi2}pi+ ${D2_pi3}pi- {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${Kp}K+]CC) || (${Bp}[B+ => ${D1}(D~0 ==> ${D1_pi1}pi- ${D1_pi2}pi+ ${D1_pi3}pi- {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${D2}(D*(2007)0 ==> ${D2_pi1}pi+ ${D2_pi2}pi- ${D2_pi3}pi+ {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${Kp}K+]CC)')
mcdtt_BuD0starD0starKp.setDescriptorTemplate('(${Bp}[B+ => ${D1}(D*(2007)~0 ==> ${D1_pi1}pi+ ${D1_pi2}pi- ${D1_pi3}pi+ {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${D2}(D*(2007)0 ==> ${D2_pi1}pi- ${D2_pi2}pi+ ${D2_pi3}pi- {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${Kp}K+]CC) || (${Bp}[B+ => ${D1}(D*(2007)~0 ==> ${D1_pi1}pi- ${D1_pi2}pi+ ${D1_pi3}pi- {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${D2}(D*(2007)0 ==> ${D2_pi1}pi+ ${D2_pi2}pi- ${D2_pi3}pi+ {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${Kp}K+]CC)')

# B+ -> D+D- K+
mcdtt_BuDpDmKp.setDescriptorTemplate('(${Bp}[B+ => ${D1}(D+ ==> ${D1_pi1}pi+ ${D1_pi2}pi- ${D1_pi3}pi+ {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${D2}(D- ==> ${D2_pi1}pi- ${D2_pi2}pi+ ${D2_pi3}pi- {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${Kp}K+]CC) || (${Bp}[B+ => ${D1}(D+ ==> ${D1_pi1}pi- ${D1_pi2}pi+ ${D1_pi3}pi- {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${D2}(D- ==> ${D2_pi1}pi+ ${D2_pi2}pi- ${D2_pi3}pi+ {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${Kp}K+]CC)')
mcdtt_BuDpstarDmKp.setDescriptorTemplate('(${Bp}[B+ => ${D1}(D*(2010)+ ==> ${D1_pi1}pi+ ${D1_pi2}pi- ${D1_pi3}pi+ {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${D2}(D- ==> ${D2_pi1}pi- ${D2_pi2}pi+ ${D2_pi3}pi- {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${Kp}K+]CC) || (${Bp}[B+ => ${D1}(D*(2010)+ ==> ${D1_pi1}pi- ${D1_pi2}pi+ ${D1_pi3}pi- {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${D2}(D- ==> ${D2_pi1}pi+ ${D2_pi2}pi- ${D2_pi3}pi+ {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${Kp}K+]CC)')
mcdtt_BuDpDmstarKp.setDescriptorTemplate('(${Bp}[B+ => ${D1}(D+ ==> ${D1_pi1}pi+ ${D1_pi2}pi- ${D1_pi3}pi+ {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${D2}(D*(2010)- ==> ${D2_pi1}pi- ${D2_pi2}pi+ ${D2_pi3}pi- {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${Kp}K+]CC) || (${Bp}[B+ => ${D1}(D+ ==> ${D1_pi1}pi- ${D1_pi2}pi+ ${D1_pi3}pi- {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${D2}(D*(2010)- ==> ${D2_pi1}pi+ ${D2_pi2}pi- ${D2_pi3}pi+ {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${Kp}K+]CC)')
mcdtt_BuDpstarDmstarKp.setDescriptorTemplate('(${Bp}[B+ => ${D1}(D*(2010)+ ==> ${D1_pi1}pi+ ${D1_pi2}pi- ${D1_pi3}pi+ {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${D2}(D*(2010)- ==> ${D2_pi1}pi- ${D2_pi2}pi+ ${D2_pi3}pi- {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${Kp}K+]CC) || (${Bp}[B+ => ${D1}(D*(2010)+ ==> ${D1_pi1}pi- ${D1_pi2}pi+ ${D1_pi3}pi- {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${D2}(D*(2010)- ==> ${D2_pi1}pi+ ${D2_pi2}pi- ${D2_pi3}pi+ {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${Kp}K+]CC)')

# B+ -> Ds+Ds- K+
mcdtt_BuDsDsKp.setDescriptorTemplate('(${Bp}[B+ => ${D1}(D_s+ ==> ${D1_pi1}pi+ ${D1_pi2}pi- ${D1_pi3}pi+ {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${D2}(D_s- ==> ${D2_pi1}pi- ${D2_pi2}pi+ ${D2_pi3}pi- {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${Kp}K+]CC) || (${Bp}[B+ => ${D1}(D_s+ ==> ${D1_pi1}pi- ${D1_pi2}pi+ ${D1_pi3}pi- {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${D2}(D_s- ==> ${D2_pi1}pi+ ${D2_pi2}pi- ${D2_pi3}pi+ {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${Kp}K+]CC)')
mcdtt_BuDsstarDsKp.setDescriptorTemplate('(${Bp}[B+ => ${D1}(D*_s+ ==> ${D1_pi1}pi+ ${D1_pi2}pi- ${D1_pi3}pi+ {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${D2}(D_s- ==> ${D2_pi1}pi- ${D2_pi2}pi+ ${D2_pi3}pi- {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${Kp}K+]CC) || (${Bp}[B+ => ${D1}(D*_s+ ==> ${D1_pi1}pi- ${D1_pi2}pi+ ${D1_pi3}pi- {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${D2}(D_s- ==> ${D2_pi1}pi+ ${D2_pi2}pi- ${D2_pi3}pi+ {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${Kp}K+]CC)')
mcdtt_BuDsDsstarKp.setDescriptorTemplate('(${Bp}[B+ => ${D1}(D_s+ ==> ${D1_pi1}pi+ ${D1_pi2}pi- ${D1_pi3}pi+ {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${D2}(D*_s- ==> ${D2_pi1}pi- ${D2_pi2}pi+ ${D2_pi3}pi- {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${Kp}K+]CC) || (${Bp}[B+ => ${D1}(D_s+ ==> ${D1_pi1}pi- ${D1_pi2}pi+ ${D1_pi3}pi- {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${D2}(D*_s- ==> ${D2_pi1}pi+ ${D2_pi2}pi- ${D2_pi3}pi+ {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${Kp}K+]CC)')
mcdtt_BuDsstarDsstarKp.setDescriptorTemplate('(${Bp}[B+ => ${D1}(D*_s+ ==> ${D1_pi1}pi+ ${D1_pi2}pi- ${D1_pi3}pi+ {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${D2}(D*_s- ==> ${D2_pi1}pi- ${D2_pi2}pi+ ${D2_pi3}pi- {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${Kp}K+]CC) || (${Bp}[B+ => ${D1}(D*_s+ ==> ${D1_pi1}pi- ${D1_pi2}pi+ ${D1_pi3}pi- {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${D2}(D*_s- ==> ${D2_pi1}pi+ ${D2_pi2}pi- ${D2_pi3}pi+ {X} {X} {X} {X} {X} {X} {X} {X} {X} {X}) ${Kp}K+]CC)')

# TUPLE TOOLS
ToolList = ['TupleToolTrackInfo', 'TupleToolVtxIsoln', 'TupleToolPrimaries', 'TupleToolCorrectedMass', 'TupleToolPropertime', 'TupleToolAngles', 'TupleToolL0Calo', 'TupleToolRecoStats']
dtt.ToolList += ToolList
if not isMC:
	dtt_SS.ToolList += ToolList
if isMC:
	dtt.ToolList += ['TupleToolMCTruth','TupleToolMCBackgroundInfo']
	dtt.addTool(TupleToolMCTruth())
	dtt.TupleToolMCTruth.ToolList += ['MCTupleToolHierarchy']

# Adds several tupletools, including DTF
addTools(dtt, stream, line, year)
if not isMC:
	addTools(dtt_SS, stream, line_SS, year)

# from Configurables import PrintMCTree, PrintMCDecayTreeTool
# mctree = PrintMCTree("PrintTrueBp")
# mctree.addTool( PrintMCDecayTreeTool, name = "PrintMC" )
# mctree.PrintMC.Information = "Name M P Px Py Pz Pt"
# mctree.ParticleNames = [ "B+", "B-" ]

# from Configurables import PrintDecayTree, PrintDecayTreeTool
# tree = PrintDecayTree("tree")
# tree.Inputs = ['/Phys/{0}/Particles'.format(line)]
# tree.addTool(PrintDecayTreeTool, name="PrintDecay")
# tree.PrintDecay.Information = "Name M P"

GENToolList = ['MCTupleToolKinematic', 'MCTupleToolHierarchy', 'TupleToolEventInfo']

if isMC:
	mcdtt.ToolList += GENToolList

	mcdtt_BuD0D0Kp.ToolList += GENToolList
	mcdtt_BuD0starD0Kp.ToolList += GENToolList
	mcdtt_BuD0D0starKp.ToolList += GENToolList
	mcdtt_BuD0starD0starKp.ToolList += GENToolList

	mcdtt_BuDpDmKp.ToolList += GENToolList
	mcdtt_BuDpstarDmKp.ToolList += GENToolList
	mcdtt_BuDpDmstarKp.ToolList += GENToolList
	mcdtt_BuDpstarDmstarKp.ToolList += GENToolList

	mcdtt_BuDsDsKp.ToolList += GENToolList
	mcdtt_BuDsstarDsKp.ToolList += GENToolList
	mcdtt_BuDsDsstarKp.ToolList += GENToolList
	mcdtt_BuDsstarDsstarKp.ToolList += GENToolList

# Stripping filter
from PhysConf.Filters import LoKi_Filters
fltrs = LoKi_Filters (
    STRIP_Code = "HLT_PASS_RE('StrippingB2KTauTau(Line|SSLine)Decision')"
	# HLT1_Code = "(HLT_PASS_RE('Hlt1TrackMVADecision')) | (HLT_PASS_RE('Hlt1TwoTrackMVADecision'))",
	# HLT2_Code = "(HLT_PASS_RE('Hlt2Topo2BodyDecision')) | (HLT_PASS_RE('Hlt2Topo3BodyDecision')) | (HLT_PASS_RE('Hlt2Topo4BodyDecision'))"
)
if not isMC:
	DaVinci().EventPreFilters = fltrs.filters('Filters')

# PV checker -> TupleToolPrimaries requires there to be a PV in the event, otherwise it throws an error
DaVinci().UserAlgorithms += [CheckPV()]

if applySel:
	DaVinci().UserAlgorithms += [b_Filter]
	if not isMC:
		DaVinci().UserAlgorithms += [b_Filter_SS]

DaVinci().UserAlgorithms += [dtt]

if isMC:
	DaVinci().UserAlgorithms += [mcdtt]
	DaVinci().UserAlgorithms += [mcdtt_BuD0D0Kp, mcdtt_BuD0starD0Kp, mcdtt_BuD0D0starKp, mcdtt_BuD0starD0starKp]
	DaVinci().UserAlgorithms += [mcdtt_BuDpDmKp, mcdtt_BuDpstarDmKp, mcdtt_BuDpDmstarKp, mcdtt_BuDpstarDmstarKp]
	DaVinci().UserAlgorithms += [mcdtt_BuDsDsKp, mcdtt_BuDsstarDsKp, mcdtt_BuDsDsstarKp, mcdtt_BuDsstarDsstarKp]
else:
    DaVinci().UserAlgorithms += [dtt_SS]

# Use the local input data
from GaudiConf import IOHelper
filename = '/panfs/felician/SimulationJobs/12693500/2016/Sim10c_ReDecay/{0}/2016_{0}.txt'.format(pol)
with open(filename) as file:
    input_files = [line.rstrip() for line in file]
IOHelper('ROOT').inputFiles(input_files, clear=True)