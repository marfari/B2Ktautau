from Configurables import DecayTreeTuple
from DecayTreeTuple.Configuration import *
from Configurables import LoKi__Hybrid__TupleTool
from Configurables import CheckPV
from Configurables import DaVinci
from PhysConf.Selections import FilterSelection
from PhysSelPython.Wrappers import AutomaticData

isMC = False
year = '2018'
pol = 'MagDown'
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
		# Vertex isolation variables
		# B+
		"VTXISONUMVTX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_VertexIsoInfo', 'VTXISONUMVTX', -100)".format(stream,line),
		"VTXISODCHI2ONETRACK_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_VertexIsoInfo', 'VTXISODCHI2ONETRACK', -100)".format(stream,line),
		"VTXISODCHI2MASSONETRACK_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_VertexIsoInfo', 'VTXISODCHI2MASSONETRACK', -100)".format(stream,line),
		"VTXISODCHI2TWOTRACK_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_VertexIsoInfo', 'VTXISODCHI2TWOTRACK', -100)".format(stream,line),
		"VTXISODCHI2MASSTWOTRACK_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_VertexIsoInfo', 'VTXISODCHI2MASSTWOTRACK', -100)".format(stream,line),

		# D+
		"VTXISONUMVTX_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_VertexIsoInfo', 'VTXISONUMVTX', -100)".format(stream,line),
		"VTXISODCHI2ONETRACK_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_VertexIsoInfo', 'VTXISODCHI2ONETRACK', -100)".format(stream,line),
		"VTXISODCHI2MASSONETRACK_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_VertexIsoInfo', 'VTXISODCHI2MASSONETRACK', -100)".format(stream,line),
		"VTXISODCHI2TWOTRACK_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_VertexIsoInfo', 'VTXISODCHI2TWOTRACK', -100)".format(stream,line),
		"VTXISODCHI2MASSTWOTRACK_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_VertexIsoInfo', 'VTXISODCHI2MASSTWOTRACK', -100)".format(stream,line),

		# D-
		"VTXISONUMVTX_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_VertexIsoInfo', 'VTXISONUMVTX', -100)".format(stream,line),
		"VTXISODCHI2ONETRACK_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_VertexIsoInfo', 'VTXISODCHI2ONETRACK', -100)".format(stream,line),
		"VTXISODCHI2MASSONETRACK_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_VertexIsoInfo', 'VTXISODCHI2MASSONETRACK', -100)".format(stream,line),
		"VTXISODCHI2TWOTRACK_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_VertexIsoInfo', 'VTXISODCHI2TWOTRACK', -100)".format(stream,line),
		"VTXISODCHI2MASSTWOTRACK_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_VertexIsoInfo', 'VTXISODCHI2MASSTWOTRACK', -100)".format(stream,line),

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
		"CC_05_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone05', 'CC_PX', -100.)".format(stream, line),
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
		"CC_05_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone05', 'CC_PX', -100.)".format(stream, line),
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

		# D+ (cone size = 0.5)
		"CC_05_ANGLE_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'CC_ANGLE', -100.)".format(stream, line),
		"CC_05_MULT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'CC_MULT', -100.)".format(stream, line),
		"CC_05_SPT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'CC_SPT', -100.)".format(stream, line),
		"CC_05_VPT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'CC_VPT', -100.)".format(stream, line),
		"CC_05_PX_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'CC_PX', -100.)".format(stream, line),
		"CC_05_PY_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'CC_PY', -100.)".format(stream, line),
		"CC_05_PZ_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'CC_PZ', -100.)".format(stream, line),
		"CC_05_PASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'CC_PASYM', -100.)".format(stream, line),
		"CC_05_PTASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'CC_PTASYM', -100.)".format(stream, line),
		"CC_05_PXASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'CC_PXASYM', -100.)".format(stream, line),
		"CC_05_PYASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'CC_PYASYM', -100.)".format(stream, line),
		"CC_05_PZASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'CC_PZASYM', -100.)".format(stream, line),
		"CC_05_DELTAETA_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'CC_DELTAETA', -100.)".format(stream, line),
		"CC_05_DELTAPHI_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'CC_DELTAPHI', -100.)".format(stream, line),
		"CC_05_PX_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'CC_PX', -100.)".format(stream, line),
		"CC_05_IT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'CC_IT', -100.)".format(stream, line),
		"CC_05_MAXPT_Q_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'CC_MAXPT_Q', -100.)".format(stream, line),
		"CC_05_MAXPT_PT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'CC_MAXPT_PT', -100.)".format(stream, line),
		"CC_05_MAXPT_PX_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'CC_MAXPT_PX', -100.)".format(stream, line),
		"CC_05_MAXPT_PY_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'CC_MAXPT_PY', -100.)".format(stream, line),
		"CC_05_MAXPT_PZ_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'CC_MAXPT_PZ', -100.)".format(stream, line),
		"CC_05_MAXPT_PE_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'CC_MAXPT_PE', -100.)".format(stream, line),
		"NC_05_ANGLE_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'NC_ANGLE', -100.)".format(stream, line),
		"NC_05_MULT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'NC_MULT', -100.)".format(stream, line),
		"NC_05_SPT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'NC_SPT', -100.)".format(stream, line),
		"NC_05_VPT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'NC_VPT', -100.)".format(stream, line),
		"NC_05_PX_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'NC_PX', -100.)".format(stream, line),
		"NC_05_PY_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'NC_PY', -100.)".format(stream, line),
		"NC_05_PZ_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'NC_PZ', -100.)".format(stream, line),
		"NC_05_PASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'NC_PASYM', -100.)".format(stream, line),
		"NC_05_PTASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'NC_PTASYM', -100.)".format(stream, line),
		"NC_05_PXASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'NC_PXASYM', -100.)".format(stream, line),
		"NC_05_PYASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'NC_PYASYM', -100.)".format(stream, line),
		"NC_05_PZASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'NC_PZASYM', -100.)".format(stream, line),
		"NC_05_DELTAETA_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'NC_DELTAETA', -100.)".format(stream, line),
		"NC_05_DELTAPHI_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'NC_DELTAPHI', -100.)".format(stream, line),
		"NC_05_IT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'NC_IT', -100.)".format(stream, line),
		"NC_05_MAXPT_PT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'NC_MAXPT_PT', -100.)".format(stream, line),
		"NC_05_MAXPT_PX_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'NC_MAXPT_PX', -100.)".format(stream, line),
		"NC_05_MAXPT_PY_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'NC_MAXPT_PY', -100.)".format(stream, line),
		"NC_05_MAXPT_PZ_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'NC_MAXPT_PZ', -100.)".format(stream, line),
		"CCNC_05_IT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone05', 'CCNC_IT', -100.)".format(stream, line),

		# D- (cone size = 0.5)
		"CC_05_ANGLE_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'CC_ANGLE', -100.)".format(stream, line),
		"CC_05_MULT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'CC_MULT', -100.)".format(stream, line),
		"CC_05_SPT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'CC_SPT', -100.)".format(stream, line),
		"CC_05_VPT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'CC_VPT', -100.)".format(stream, line),
		"CC_05_PX_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'CC_PX', -100.)".format(stream, line),
		"CC_05_PY_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'CC_PY', -100.)".format(stream, line),
		"CC_05_PZ_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'CC_PZ', -100.)".format(stream, line),
		"CC_05_PASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'CC_PASYM', -100.)".format(stream, line),
		"CC_05_PTASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'CC_PTASYM', -100.)".format(stream, line),
		"CC_05_PXASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'CC_PXASYM', -100.)".format(stream, line),
		"CC_05_PYASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'CC_PYASYM', -100.)".format(stream, line),
		"CC_05_PZASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'CC_PZASYM', -100.)".format(stream, line),
		"CC_05_DELTAETA_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'CC_DELTAETA', -100.)".format(stream, line),
		"CC_05_DELTAPHI_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'CC_DELTAPHI', -100.)".format(stream, line),
		"CC_05_PX_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'CC_PX', -100.)".format(stream, line),
		"CC_05_IT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'CC_IT', -100.)".format(stream, line),
		"CC_05_MAXPT_Q_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'CC_MAXPT_Q', -100.)".format(stream, line),
		"CC_05_MAXPT_PT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'CC_MAXPT_PT', -100.)".format(stream, line),
		"CC_05_MAXPT_PX_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'CC_MAXPT_PX', -100.)".format(stream, line),
		"CC_05_MAXPT_PY_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'CC_MAXPT_PY', -100.)".format(stream, line),
		"CC_05_MAXPT_PZ_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'CC_MAXPT_PZ', -100.)".format(stream, line),
		"CC_05_MAXPT_PE_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'CC_MAXPT_PE', -100.)".format(stream, line),
		"NC_05_ANGLE_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'NC_ANGLE', -100.)".format(stream, line),
		"NC_05_MULT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'NC_MULT', -100.)".format(stream, line),
		"NC_05_SPT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'NC_SPT', -100.)".format(stream, line),
		"NC_05_VPT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'NC_VPT', -100.)".format(stream, line),
		"NC_05_PX_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'NC_PX', -100.)".format(stream, line),
		"NC_05_PY_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'NC_PY', -100.)".format(stream, line),
		"NC_05_PZ_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'NC_PZ', -100.)".format(stream, line),
		"NC_05_PASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'NC_PASYM', -100.)".format(stream, line),
		"NC_05_PTASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'NC_PTASYM', -100.)".format(stream, line),
		"NC_05_PXASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'NC_PXASYM', -100.)".format(stream, line),
		"NC_05_PYASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'NC_PYASYM', -100.)".format(stream, line),
		"NC_05_PZASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'NC_PZASYM', -100.)".format(stream, line),
		"NC_05_DELTAETA_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'NC_DELTAETA', -100.)".format(stream, line),
		"NC_05_DELTAPHI_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'NC_DELTAPHI', -100.)".format(stream, line),
		"NC_05_IT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'NC_IT', -100.)".format(stream, line),
		"NC_05_MAXPT_PT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'NC_MAXPT_PT', -100.)".format(stream, line),
		"NC_05_MAXPT_PX_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'NC_MAXPT_PX', -100.)".format(stream, line),
		"NC_05_MAXPT_PY_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'NC_MAXPT_PY', -100.)".format(stream, line),
		"NC_05_MAXPT_PZ_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'NC_MAXPT_PZ', -100.)".format(stream, line),
		"CCNC_05_IT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone05', 'CCNC_IT', -100.)".format(stream, line),

		# B+ (cone size = 1.0)
		"CC_10_ANGLE_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_ANGLE', -100.)".format(stream, line),
		"CC_10_MULT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_MULT', -100.)".format(stream, line),
		"CC_10_SPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_SPT', -100.)".format(stream, line),
		"CC_10_VPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_VPT', -100.)".format(stream, line),
		"CC_10_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_PX', -100.)".format(stream, line),
		"CC_10_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_PY', -100.)".format(stream, line),
		"CC_10_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_PZ', -100.)".format(stream, line),
		"CC_10_PASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_PASYM', -100.)".format(stream, line),
		"CC_10_PTASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_PTASYM', -100.)".format(stream, line),
		"CC_10_PXASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_PXASYM', -100.)".format(stream, line),
		"CC_10_PYASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_PYASYM', -100.)".format(stream, line),
		"CC_10_PZASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_PZASYM', -100.)".format(stream, line),
		"CC_10_DELTAETA_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_DELTAETA', -100.)".format(stream, line),
		"CC_10_DELTAPHI_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_DELTAPHI', -100.)".format(stream, line),
		"CC_10_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_PX', -100.)".format(stream, line),
		"CC_10_IT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_IT', -100.)".format(stream, line),
		"CC_10_MAXPT_Q_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_MAXPT_Q', -100.)".format(stream, line),
		"CC_10_MAXPT_PT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_MAXPT_PT', -100.)".format(stream, line),
		"CC_10_MAXPT_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_MAXPT_PX', -100.)".format(stream, line),
		"CC_10_MAXPT_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_MAXPT_PY', -100.)".format(stream, line),
		"CC_10_MAXPT_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_MAXPT_PZ', -100.)".format(stream, line),
		"CC_10_MAXPT_PE_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CC_MAXPT_PE', -100.)".format(stream, line),
		"NC_10_ANGLE_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_ANGLE', -100.)".format(stream, line),
		"NC_10_MULT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_MULT', -100.)".format(stream, line),
		"NC_10_SPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_SPT', -100.)".format(stream, line),
		"NC_10_VPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_VPT', -100.)".format(stream, line),
		"NC_10_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_PX', -100.)".format(stream, line),
		"NC_10_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_PY', -100.)".format(stream, line),
		"NC_10_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_PZ', -100.)".format(stream, line),
		"NC_10_PASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_PASYM', -100.)".format(stream, line),
		"NC_10_PTASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_PTASYM', -100.)".format(stream, line),
		"NC_10_PXASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_PXASYM', -100.)".format(stream, line),
		"NC_10_PYASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_PYASYM', -100.)".format(stream, line),
		"NC_10_PZASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_PZASYM', -100.)".format(stream, line),
		"NC_10_DELTAETA_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_DELTAETA', -100.)".format(stream, line),
		"NC_10_DELTAPHI_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_DELTAPHI', -100.)".format(stream, line),
		"NC_10_IT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_IT', -100.)".format(stream, line),
		"NC_10_MAXPT_PT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_MAXPT_PT', -100.)".format(stream, line),
		"NC_10_MAXPT_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_MAXPT_PX', -100.)".format(stream, line),
		"NC_10_MAXPT_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_MAXPT_PY', -100.)".format(stream, line),
		"NC_10_MAXPT_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'NC_MAXPT_PZ', -100.)".format(stream, line),
		"CCNC_10_IT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone10', 'CCNC_IT', -100.)".format(stream, line),

		# K+ (cone size = 1.0)
		"CC_10_ANGLE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_ANGLE', -100.)".format(stream, line),
		"CC_10_MULT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_MULT', -100.)".format(stream, line),
		"CC_10_SPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_SPT', -100.)".format(stream, line),
		"CC_10_VPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_VPT', -100.)".format(stream, line),
		"CC_10_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_PX', -100.)".format(stream, line),
		"CC_10_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_PY', -100.)".format(stream, line),
		"CC_10_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_PZ', -100.)".format(stream, line),
		"CC_10_PASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_PASYM', -100.)".format(stream, line),
		"CC_10_PTASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_PTASYM', -100.)".format(stream, line),
		"CC_10_PXASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_PXASYM', -100.)".format(stream, line),
		"CC_10_PYASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_PYASYM', -100.)".format(stream, line),
		"CC_10_PZASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_PZASYM', -100.)".format(stream, line),
		"CC_10_DELTAETA_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_DELTAETA', -100.)".format(stream, line),
		"CC_10_DELTAPHI_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_DELTAPHI', -100.)".format(stream, line),
		"CC_10_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_PX', -100.)".format(stream, line),
		"CC_10_IT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_IT', -100.)".format(stream, line),
		"CC_10_MAXPT_Q_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_MAXPT_Q', -100.)".format(stream, line),
		"CC_10_MAXPT_PT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_MAXPT_PT', -100.)".format(stream, line),
		"CC_10_MAXPT_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_MAXPT_PX', -100.)".format(stream, line),
		"CC_10_MAXPT_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_MAXPT_PY', -100.)".format(stream, line),
		"CC_10_MAXPT_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_MAXPT_PZ', -100.)".format(stream, line),
		"CC_10_MAXPT_PE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CC_MAXPT_PE', -100.)".format(stream, line),
		"NC_10_ANGLE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_ANGLE', -100.)".format(stream, line),
		"NC_10_MULT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_MULT', -100.)".format(stream, line),
		"NC_10_SPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_SPT', -100.)".format(stream, line),
		"NC_10_VPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_VPT', -100.)".format(stream, line),
		"NC_10_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_PX', -100.)".format(stream, line),
		"NC_10_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_PY', -100.)".format(stream, line),
		"NC_10_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_PZ', -100.)".format(stream, line),
		"NC_10_PASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_PASYM', -100.)".format(stream, line),
		"NC_10_PTASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_PTASYM', -100.)".format(stream, line),
		"NC_10_PXASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_PXASYM', -100.)".format(stream, line),
		"NC_10_PYASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_PYASYM', -100.)".format(stream, line),
		"NC_10_PZASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_PZASYM', -100.)".format(stream, line),
		"NC_10_DELTAETA_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_DELTAETA', -100.)".format(stream, line),
		"NC_10_DELTAPHI_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_DELTAPHI', -100.)".format(stream, line),
		"NC_10_IT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_IT', -100.)".format(stream, line),
		"NC_10_MAXPT_PT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_MAXPT_PT', -100.)".format(stream, line),
		"NC_10_MAXPT_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_MAXPT_PX', -100.)".format(stream, line),
		"NC_10_MAXPT_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_MAXPT_PY', -100.)".format(stream, line),
		"NC_10_MAXPT_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'NC_MAXPT_PZ', -100.)".format(stream, line),
		"CCNC_10_IT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone10', 'CCNC_IT', -100.)".format(stream, line),

		# D+ (cone size = 1.0)
		"CC_10_ANGLE_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'CC_ANGLE', -100.)".format(stream, line),
		"CC_10_MULT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'CC_MULT', -100.)".format(stream, line),
		"CC_10_SPT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'CC_SPT', -100.)".format(stream, line),
		"CC_10_VPT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'CC_VPT', -100.)".format(stream, line),
		"CC_10_PX_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'CC_PX', -100.)".format(stream, line),
		"CC_10_PY_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'CC_PY', -100.)".format(stream, line),
		"CC_10_PZ_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'CC_PZ', -100.)".format(stream, line),
		"CC_10_PASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'CC_PASYM', -100.)".format(stream, line),
		"CC_10_PTASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'CC_PTASYM', -100.)".format(stream, line),
		"CC_10_PXASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'CC_PXASYM', -100.)".format(stream, line),
		"CC_10_PYASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'CC_PYASYM', -100.)".format(stream, line),
		"CC_10_PZASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'CC_PZASYM', -100.)".format(stream, line),
		"CC_10_DELTAETA_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'CC_DELTAETA', -100.)".format(stream, line),
		"CC_10_DELTAPHI_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'CC_DELTAPHI', -100.)".format(stream, line),
		"CC_10_PX_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'CC_PX', -100.)".format(stream, line),
		"CC_10_IT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'CC_IT', -100.)".format(stream, line),
		"CC_10_MAXPT_Q_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'CC_MAXPT_Q', -100.)".format(stream, line),
		"CC_10_MAXPT_PT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'CC_MAXPT_PT', -100.)".format(stream, line),
		"CC_10_MAXPT_PX_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'CC_MAXPT_PX', -100.)".format(stream, line),
		"CC_10_MAXPT_PY_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'CC_MAXPT_PY', -100.)".format(stream, line),
		"CC_10_MAXPT_PZ_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'CC_MAXPT_PZ', -100.)".format(stream, line),
		"CC_10_MAXPT_PE_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'CC_MAXPT_PE', -100.)".format(stream, line),
		"NC_10_ANGLE_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'NC_ANGLE', -100.)".format(stream, line),
		"NC_10_MULT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'NC_MULT', -100.)".format(stream, line),
		"NC_10_SPT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'NC_SPT', -100.)".format(stream, line),
		"NC_10_VPT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'NC_VPT', -100.)".format(stream, line),
		"NC_10_PX_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'NC_PX', -100.)".format(stream, line),
		"NC_10_PY_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'NC_PY', -100.)".format(stream, line),
		"NC_10_PZ_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'NC_PZ', -100.)".format(stream, line),
		"NC_10_PASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'NC_PASYM', -100.)".format(stream, line),
		"NC_10_PTASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'NC_PTASYM', -100.)".format(stream, line),
		"NC_10_PXASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'NC_PXASYM', -100.)".format(stream, line),
		"NC_10_PYASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'NC_PYASYM', -100.)".format(stream, line),
		"NC_10_PZASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'NC_PZASYM', -100.)".format(stream, line),
		"NC_10_DELTAETA_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'NC_DELTAETA', -100.)".format(stream, line),
		"NC_10_DELTAPHI_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'NC_DELTAPHI', -100.)".format(stream, line),
		"NC_10_IT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'NC_IT', -100.)".format(stream, line),
		"NC_10_MAXPT_PT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'NC_MAXPT_PT', -100.)".format(stream, line),
		"NC_10_MAXPT_PX_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'NC_MAXPT_PX', -100.)".format(stream, line),
		"NC_10_MAXPT_PY_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'NC_MAXPT_PY', -100.)".format(stream, line),
		"NC_10_MAXPT_PZ_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'NC_MAXPT_PZ', -100.)".format(stream, line),
		"CCNC_10_IT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone10', 'CCNC_IT', -100.)".format(stream, line),

		# D- (cone size = 1.0)
		"CC_10_ANGLE_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'CC_ANGLE', -100.)".format(stream, line),
		"CC_10_MULT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'CC_MULT', -100.)".format(stream, line),
		"CC_10_SPT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'CC_SPT', -100.)".format(stream, line),
		"CC_10_VPT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'CC_VPT', -100.)".format(stream, line),
		"CC_10_PX_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'CC_PX', -100.)".format(stream, line),
		"CC_10_PY_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'CC_PY', -100.)".format(stream, line),
		"CC_10_PZ_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'CC_PZ', -100.)".format(stream, line),
		"CC_10_PASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'CC_PASYM', -100.)".format(stream, line),
		"CC_10_PTASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'CC_PTASYM', -100.)".format(stream, line),
		"CC_10_PXASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'CC_PXASYM', -100.)".format(stream, line),
		"CC_10_PYASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'CC_PYASYM', -100.)".format(stream, line),
		"CC_10_PZASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'CC_PZASYM', -100.)".format(stream, line),
		"CC_10_DELTAETA_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'CC_DELTAETA', -100.)".format(stream, line),
		"CC_10_DELTAPHI_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'CC_DELTAPHI', -100.)".format(stream, line),
		"CC_10_PX_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'CC_PX', -100.)".format(stream, line),
		"CC_10_IT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'CC_IT', -100.)".format(stream, line),
		"CC_10_MAXPT_Q_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'CC_MAXPT_Q', -100.)".format(stream, line),
		"CC_10_MAXPT_PT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'CC_MAXPT_PT', -100.)".format(stream, line),
		"CC_10_MAXPT_PX_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'CC_MAXPT_PX', -100.)".format(stream, line),
		"CC_10_MAXPT_PY_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'CC_MAXPT_PY', -100.)".format(stream, line),
		"CC_10_MAXPT_PZ_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'CC_MAXPT_PZ', -100.)".format(stream, line),
		"CC_10_MAXPT_PE_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'CC_MAXPT_PE', -100.)".format(stream, line),
		"NC_10_ANGLE_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'NC_ANGLE', -100.)".format(stream, line),
		"NC_10_MULT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'NC_MULT', -100.)".format(stream, line),
		"NC_10_SPT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'NC_SPT', -100.)".format(stream, line),
		"NC_10_VPT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'NC_VPT', -100.)".format(stream, line),
		"NC_10_PX_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'NC_PX', -100.)".format(stream, line),
		"NC_10_PY_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'NC_PY', -100.)".format(stream, line),
		"NC_10_PZ_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'NC_PZ', -100.)".format(stream, line),
		"NC_10_PASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'NC_PASYM', -100.)".format(stream, line),
		"NC_10_PTASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'NC_PTASYM', -100.)".format(stream, line),
		"NC_10_PXASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'NC_PXASYM', -100.)".format(stream, line),
		"NC_10_PYASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'NC_PYASYM', -100.)".format(stream, line),
		"NC_10_PZASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'NC_PZASYM', -100.)".format(stream, line),
		"NC_10_DELTAETA_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'NC_DELTAETA', -100.)".format(stream, line),
		"NC_10_DELTAPHI_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'NC_DELTAPHI', -100.)".format(stream, line),
		"NC_10_IT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'NC_IT', -100.)".format(stream, line),
		"NC_10_MAXPT_PT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'NC_MAXPT_PT', -100.)".format(stream, line),
		"NC_10_MAXPT_PX_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'NC_MAXPT_PX', -100.)".format(stream, line),
		"NC_10_MAXPT_PY_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'NC_MAXPT_PY', -100.)".format(stream, line),
		"NC_10_MAXPT_PZ_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'NC_MAXPT_PZ', -100.)".format(stream, line),
		"CCNC_10_IT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone10', 'CCNC_IT', -100.)".format(stream, line),

		# B+ (cone size = 1.5)
		"CC_15_ANGLE_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_ANGLE', -100.)".format(stream, line),
		"CC_15_MULT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_MULT', -100.)".format(stream, line),
		"CC_15_SPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_SPT', -100.)".format(stream, line),
		"CC_15_VPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_VPT', -100.)".format(stream, line),
		"CC_15_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_PX', -100.)".format(stream, line),
		"CC_15_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_PY', -100.)".format(stream, line),
		"CC_15_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_PZ', -100.)".format(stream, line),
		"CC_15_PASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_PASYM', -100.)".format(stream, line),
		"CC_15_PTASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_PTASYM', -100.)".format(stream, line),
		"CC_15_PXASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_PXASYM', -100.)".format(stream, line),
		"CC_15_PYASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_PYASYM', -100.)".format(stream, line),
		"CC_15_PZASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_PZASYM', -100.)".format(stream, line),
		"CC_15_DELTAETA_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_DELTAETA', -100.)".format(stream, line),
		"CC_15_DELTAPHI_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_DELTAPHI', -100.)".format(stream, line),
		"CC_15_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_PX', -100.)".format(stream, line),
		"CC_15_IT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_IT', -100.)".format(stream, line),
		"CC_15_MAXPT_Q_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_MAXPT_Q', -100.)".format(stream, line),
		"CC_15_MAXPT_PT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_MAXPT_PT', -100.)".format(stream, line),
		"CC_15_MAXPT_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_MAXPT_PX', -100.)".format(stream, line),
		"CC_15_MAXPT_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_MAXPT_PY', -100.)".format(stream, line),
		"CC_15_MAXPT_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_MAXPT_PZ', -100.)".format(stream, line),
		"CC_15_MAXPT_PE_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CC_MAXPT_PE', -100.)".format(stream, line),
		"NC_15_ANGLE_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_ANGLE', -100.)".format(stream, line),
		"NC_15_MULT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_MULT', -100.)".format(stream, line),
		"NC_15_SPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_SPT', -100.)".format(stream, line),
		"NC_15_VPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_VPT', -100.)".format(stream, line),
		"NC_15_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_PX', -100.)".format(stream, line),
		"NC_15_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_PY', -100.)".format(stream, line),
		"NC_15_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_PZ', -100.)".format(stream, line),
		"NC_15_PASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_PASYM', -100.)".format(stream, line),
		"NC_15_PTASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_PTASYM', -100.)".format(stream, line),
		"NC_15_PXASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_PXASYM', -100.)".format(stream, line),
		"NC_15_PYASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_PYASYM', -100.)".format(stream, line),
		"NC_15_PZASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_PZASYM', -100.)".format(stream, line),
		"NC_15_DELTAETA_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_DELTAETA', -100.)".format(stream, line),
		"NC_15_DELTAPHI_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_DELTAPHI', -100.)".format(stream, line),
		"NC_15_IT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_IT', -100.)".format(stream, line),
		"NC_15_MAXPT_PT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_MAXPT_PT', -100.)".format(stream, line),
		"NC_15_MAXPT_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_MAXPT_PX', -100.)".format(stream, line),
		"NC_15_MAXPT_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_MAXPT_PY', -100.)".format(stream, line),
		"NC_15_MAXPT_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'NC_MAXPT_PZ', -100.)".format(stream, line),
		"CCNC_15_IT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone15', 'CCNC_IT', -100.)".format(stream, line),

		# K+ (cone size = 1.5)
		"CC_15_ANGLE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_ANGLE', -100.)".format(stream, line),
		"CC_15_MULT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_MULT', -100.)".format(stream, line),
		"CC_15_SPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_SPT', -100.)".format(stream, line),
		"CC_15_VPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_VPT', -100.)".format(stream, line),
		"CC_15_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_PX', -100.)".format(stream, line),
		"CC_15_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_PY', -100.)".format(stream, line),
		"CC_15_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_PZ', -100.)".format(stream, line),
		"CC_15_PASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_PASYM', -100.)".format(stream, line),
		"CC_15_PTASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_PTASYM', -100.)".format(stream, line),
		"CC_15_PXASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_PXASYM', -100.)".format(stream, line),
		"CC_15_PYASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_PYASYM', -100.)".format(stream, line),
		"CC_15_PZASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_PZASYM', -100.)".format(stream, line),
		"CC_15_DELTAETA_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_DELTAETA', -100.)".format(stream, line),
		"CC_15_DELTAPHI_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_DELTAPHI', -100.)".format(stream, line),
		"CC_15_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_PX', -100.)".format(stream, line),
		"CC_15_IT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_IT', -100.)".format(stream, line),
		"CC_15_MAXPT_Q_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_MAXPT_Q', -100.)".format(stream, line),
		"CC_15_MAXPT_PT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_MAXPT_PT', -100.)".format(stream, line),
		"CC_15_MAXPT_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_MAXPT_PX', -100.)".format(stream, line),
		"CC_15_MAXPT_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_MAXPT_PY', -100.)".format(stream, line),
		"CC_15_MAXPT_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_MAXPT_PZ', -100.)".format(stream, line),
		"CC_15_MAXPT_PE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CC_MAXPT_PE', -100.)".format(stream, line),
		"NC_15_ANGLE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_ANGLE', -100.)".format(stream, line),
		"NC_15_MULT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_MULT', -100.)".format(stream, line),
		"NC_15_SPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_SPT', -100.)".format(stream, line),
		"NC_15_VPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_VPT', -100.)".format(stream, line),
		"NC_15_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_PX', -100.)".format(stream, line),
		"NC_15_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_PY', -100.)".format(stream, line),
		"NC_15_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_PZ', -100.)".format(stream, line),
		"NC_15_PASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_PASYM', -100.)".format(stream, line),
		"NC_15_PTASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_PTASYM', -100.)".format(stream, line),
		"NC_15_PXASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_PXASYM', -100.)".format(stream, line),
		"NC_15_PYASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_PYASYM', -100.)".format(stream, line),
		"NC_15_PZASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_PZASYM', -100.)".format(stream, line),
		"NC_15_DELTAETA_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_DELTAETA', -100.)".format(stream, line),
		"NC_15_DELTAPHI_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_DELTAPHI', -100.)".format(stream, line),
		"NC_15_IT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_IT', -100.)".format(stream, line),
		"NC_15_MAXPT_PT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_MAXPT_PT', -100.)".format(stream, line),
		"NC_15_MAXPT_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_MAXPT_PX', -100.)".format(stream, line),
		"NC_15_MAXPT_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_MAXPT_PY', -100.)".format(stream, line),
		"NC_15_MAXPT_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'NC_MAXPT_PZ', -100.)".format(stream, line),
		"CCNC_15_IT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone15', 'CCNC_IT', -100.)".format(stream, line),

		# D+ (cone size = 1.5)
		"CC_15_ANGLE_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'CC_ANGLE', -100.)".format(stream, line),
		"CC_15_MULT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'CC_MULT', -100.)".format(stream, line),
		"CC_15_SPT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'CC_SPT', -100.)".format(stream, line),
		"CC_15_VPT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'CC_VPT', -100.)".format(stream, line),
		"CC_15_PX_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'CC_PX', -100.)".format(stream, line),
		"CC_15_PY_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'CC_PY', -100.)".format(stream, line),
		"CC_15_PZ_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'CC_PZ', -100.)".format(stream, line),
		"CC_15_PASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'CC_PASYM', -100.)".format(stream, line),
		"CC_15_PTASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'CC_PTASYM', -100.)".format(stream, line),
		"CC_15_PXASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'CC_PXASYM', -100.)".format(stream, line),
		"CC_15_PYASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'CC_PYASYM', -100.)".format(stream, line),
		"CC_15_PZASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'CC_PZASYM', -100.)".format(stream, line),
		"CC_15_DELTAETA_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'CC_DELTAETA', -100.)".format(stream, line),
		"CC_15_DELTAPHI_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'CC_DELTAPHI', -100.)".format(stream, line),
		"CC_15_PX_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'CC_PX', -100.)".format(stream, line),
		"CC_15_IT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'CC_IT', -100.)".format(stream, line),
		"CC_15_MAXPT_Q_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'CC_MAXPT_Q', -100.)".format(stream, line),
		"CC_15_MAXPT_PT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'CC_MAXPT_PT', -100.)".format(stream, line),
		"CC_15_MAXPT_PX_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'CC_MAXPT_PX', -100.)".format(stream, line),
		"CC_15_MAXPT_PY_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'CC_MAXPT_PY', -100.)".format(stream, line),
		"CC_15_MAXPT_PZ_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'CC_MAXPT_PZ', -100.)".format(stream, line),
		"CC_15_MAXPT_PE_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'CC_MAXPT_PE', -100.)".format(stream, line),
		"NC_15_ANGLE_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'NC_ANGLE', -100.)".format(stream, line),
		"NC_15_MULT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'NC_MULT', -100.)".format(stream, line),
		"NC_15_SPT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'NC_SPT', -100.)".format(stream, line),
		"NC_15_VPT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'NC_VPT', -100.)".format(stream, line),
		"NC_15_PX_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'NC_PX', -100.)".format(stream, line),
		"NC_15_PY_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'NC_PY', -100.)".format(stream, line),
		"NC_15_PZ_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'NC_PZ', -100.)".format(stream, line),
		"NC_15_PASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'NC_PASYM', -100.)".format(stream, line),
		"NC_15_PTASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'NC_PTASYM', -100.)".format(stream, line),
		"NC_15_PXASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'NC_PXASYM', -100.)".format(stream, line),
		"NC_15_PYASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'NC_PYASYM', -100.)".format(stream, line),
		"NC_15_PZASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'NC_PZASYM', -100.)".format(stream, line),
		"NC_15_DELTAETA_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'NC_DELTAETA', -100.)".format(stream, line),
		"NC_15_DELTAPHI_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'NC_DELTAPHI', -100.)".format(stream, line),
		"NC_15_IT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'NC_IT', -100.)".format(stream, line),
		"NC_15_MAXPT_PT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'NC_MAXPT_PT', -100.)".format(stream, line),
		"NC_15_MAXPT_PX_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'NC_MAXPT_PX', -100.)".format(stream, line),
		"NC_15_MAXPT_PY_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'NC_MAXPT_PY', -100.)".format(stream, line),
		"NC_15_MAXPT_PZ_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'NC_MAXPT_PZ', -100.)".format(stream, line),
		"CCNC_15_IT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone15', 'CCNC_IT', -100.)".format(stream, line),

		# tau- (cone size = 1.5)
		"CC_15_ANGLE_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'CC_ANGLE', -100.)".format(stream, line),
		"CC_15_MULT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'CC_MULT', -100.)".format(stream, line),
		"CC_15_SPT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'CC_SPT', -100.)".format(stream, line),
		"CC_15_VPT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'CC_VPT', -100.)".format(stream, line),
		"CC_15_PX_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'CC_PX', -100.)".format(stream, line),
		"CC_15_PY_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'CC_PY', -100.)".format(stream, line),
		"CC_15_PZ_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'CC_PZ', -100.)".format(stream, line),
		"CC_15_PASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'CC_PASYM', -100.)".format(stream, line),
		"CC_15_PTASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'CC_PTASYM', -100.)".format(stream, line),
		"CC_15_PXASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'CC_PXASYM', -100.)".format(stream, line),
		"CC_15_PYASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'CC_PYASYM', -100.)".format(stream, line),
		"CC_15_PZASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'CC_PZASYM', -100.)".format(stream, line),
		"CC_15_DELTAETA_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'CC_DELTAETA', -100.)".format(stream, line),
		"CC_15_DELTAPHI_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'CC_DELTAPHI', -100.)".format(stream, line),
		"CC_15_PX_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'CC_PX', -100.)".format(stream, line),
		"CC_15_IT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'CC_IT', -100.)".format(stream, line),
		"CC_15_MAXPT_Q_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'CC_MAXPT_Q', -100.)".format(stream, line),
		"CC_15_MAXPT_PT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'CC_MAXPT_PT', -100.)".format(stream, line),
		"CC_15_MAXPT_PX_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'CC_MAXPT_PX', -100.)".format(stream, line),
		"CC_15_MAXPT_PY_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'CC_MAXPT_PY', -100.)".format(stream, line),
		"CC_15_MAXPT_PZ_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'CC_MAXPT_PZ', -100.)".format(stream, line),
		"CC_15_MAXPT_PE_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'CC_MAXPT_PE', -100.)".format(stream, line),
		"NC_15_ANGLE_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'NC_ANGLE', -100.)".format(stream, line),
		"NC_15_MULT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'NC_MULT', -100.)".format(stream, line),
		"NC_15_SPT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'NC_SPT', -100.)".format(stream, line),
		"NC_15_VPT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'NC_VPT', -100.)".format(stream, line),
		"NC_15_PX_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'NC_PX', -100.)".format(stream, line),
		"NC_15_PY_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'NC_PY', -100.)".format(stream, line),
		"NC_15_PZ_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'NC_PZ', -100.)".format(stream, line),
		"NC_15_PASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'NC_PASYM', -100.)".format(stream, line),
		"NC_15_PTASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'NC_PTASYM', -100.)".format(stream, line),
		"NC_15_PXASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'NC_PXASYM', -100.)".format(stream, line),
		"NC_15_PYASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'NC_PYASYM', -100.)".format(stream, line),
		"NC_15_PZASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'NC_PZASYM', -100.)".format(stream, line),
		"NC_15_DELTAETA_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'NC_DELTAETA', -100.)".format(stream, line),
		"NC_15_DELTAPHI_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'NC_DELTAPHI', -100.)".format(stream, line),
		"NC_15_IT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'NC_IT', -100.)".format(stream, line),
		"NC_15_MAXPT_PT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'NC_MAXPT_PT', -100.)".format(stream, line),
		"NC_15_MAXPT_PX_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'NC_MAXPT_PX', -100.)".format(stream, line),
		"NC_15_MAXPT_PY_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'NC_MAXPT_PY', -100.)".format(stream, line),
		"NC_15_MAXPT_PZ_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'NC_MAXPT_PZ', -100.)".format(stream, line),
		"CCNC_15_IT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone15', 'CCNC_IT', -100.)".format(stream, line),

		# B+ (cone size = 2.0)
		"CC_20_ANGLE_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_ANGLE', -100.)".format(stream, line),
		"CC_20_MULT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_MULT', -100.)".format(stream, line),
		"CC_20_SPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_SPT', -100.)".format(stream, line),
		"CC_20_VPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_VPT', -100.)".format(stream, line),
		"CC_20_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_PX', -100.)".format(stream, line),
		"CC_20_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_PY', -100.)".format(stream, line),
		"CC_20_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_PZ', -100.)".format(stream, line),
		"CC_20_PASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_PASYM', -100.)".format(stream, line),
		"CC_20_PTASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_PTASYM', -100.)".format(stream, line),
		"CC_20_PXASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_PXASYM', -100.)".format(stream, line),
		"CC_20_PYASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_PYASYM', -100.)".format(stream, line),
		"CC_20_PZASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_PZASYM', -100.)".format(stream, line),
		"CC_20_DELTAETA_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_DELTAETA', -100.)".format(stream, line),
		"CC_20_DELTAPHI_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_DELTAPHI', -100.)".format(stream, line),
		"CC_20_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_PX', -100.)".format(stream, line),
		"CC_20_IT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_IT', -100.)".format(stream, line),
		"CC_20_MAXPT_Q_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_MAXPT_Q', -100.)".format(stream, line),
		"CC_20_MAXPT_PT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_MAXPT_PT', -100.)".format(stream, line),
		"CC_20_MAXPT_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_MAXPT_PX', -100.)".format(stream, line),
		"CC_20_MAXPT_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_MAXPT_PY', -100.)".format(stream, line),
		"CC_20_MAXPT_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_MAXPT_PZ', -100.)".format(stream, line),
		"CC_20_MAXPT_PE_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CC_MAXPT_PE', -100.)".format(stream, line),
		"NC_20_ANGLE_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_ANGLE', -100.)".format(stream, line),
		"NC_20_MULT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_MULT', -100.)".format(stream, line),
		"NC_20_SPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_SPT', -100.)".format(stream, line),
		"NC_20_VPT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_VPT', -100.)".format(stream, line),
		"NC_20_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_PX', -100.)".format(stream, line),
		"NC_20_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_PY', -100.)".format(stream, line),
		"NC_20_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_PZ', -100.)".format(stream, line),
		"NC_20_PASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_PASYM', -100.)".format(stream, line),
		"NC_20_PTASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_PTASYM', -100.)".format(stream, line),
		"NC_20_PXASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_PXASYM', -100.)".format(stream, line),
		"NC_20_PYASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_PYASYM', -100.)".format(stream, line),
		"NC_20_PZASYM_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_PZASYM', -100.)".format(stream, line),
		"NC_20_DELTAETA_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_DELTAETA', -100.)".format(stream, line),
		"NC_20_DELTAPHI_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_DELTAPHI', -100.)".format(stream, line),
		"NC_20_IT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_IT', -100.)".format(stream, line),
		"NC_20_MAXPT_PT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_MAXPT_PT', -100.)".format(stream, line),
		"NC_20_MAXPT_PX_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_MAXPT_PX', -100.)".format(stream, line),
		"NC_20_MAXPT_PY_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_MAXPT_PY', -100.)".format(stream, line),
		"NC_20_MAXPT_PZ_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'NC_MAXPT_PZ', -100.)".format(stream, line),
		"CCNC_20_IT_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_ConeIsoInfo_Cone20', 'CCNC_IT', -100.)".format(stream, line),

		# K+ (cone size = 2.0)
		"CC_20_ANGLE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_ANGLE', -100.)".format(stream, line),
		"CC_20_MULT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_MULT', -100.)".format(stream, line),
		"CC_20_SPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_SPT', -100.)".format(stream, line),
		"CC_20_VPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_VPT', -100.)".format(stream, line),
		"CC_20_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_PX', -100.)".format(stream, line),
		"CC_20_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_PY', -100.)".format(stream, line),
		"CC_20_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_PZ', -100.)".format(stream, line),
		"CC_20_PASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_PASYM', -100.)".format(stream, line),
		"CC_20_PTASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_PTASYM', -100.)".format(stream, line),
		"CC_20_PXASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_PXASYM', -100.)".format(stream, line),
		"CC_20_PYASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_PYASYM', -100.)".format(stream, line),
		"CC_20_PZASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_PZASYM', -100.)".format(stream, line),
		"CC_20_DELTAETA_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_DELTAETA', -100.)".format(stream, line),
		"CC_20_DELTAPHI_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_DELTAPHI', -100.)".format(stream, line),
		"CC_20_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_PX', -100.)".format(stream, line),
		"CC_20_IT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_IT', -100.)".format(stream, line),
		"CC_20_MAXPT_Q_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_MAXPT_Q', -100.)".format(stream, line),
		"CC_20_MAXPT_PT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_MAXPT_PT', -100.)".format(stream, line),
		"CC_20_MAXPT_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_MAXPT_PX', -100.)".format(stream, line),
		"CC_20_MAXPT_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_MAXPT_PY', -100.)".format(stream, line),
		"CC_20_MAXPT_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_MAXPT_PZ', -100.)".format(stream, line),
		"CC_20_MAXPT_PE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CC_MAXPT_PE', -100.)".format(stream, line),
		"NC_20_ANGLE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_ANGLE', -100.)".format(stream, line),
		"NC_20_MULT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_MULT', -100.)".format(stream, line),
		"NC_20_SPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_SPT', -100.)".format(stream, line),
		"NC_20_VPT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_VPT', -100.)".format(stream, line),
		"NC_20_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_PX', -100.)".format(stream, line),
		"NC_20_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_PY', -100.)".format(stream, line),
		"NC_20_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_PZ', -100.)".format(stream, line),
		"NC_20_PASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_PASYM', -100.)".format(stream, line),
		"NC_20_PTASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_PTASYM', -100.)".format(stream, line),
		"NC_20_PXASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_PXASYM', -100.)".format(stream, line),
		"NC_20_PYASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_PYASYM', -100.)".format(stream, line),
		"NC_20_PZASYM_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_PZASYM', -100.)".format(stream, line),
		"NC_20_DELTAETA_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_DELTAETA', -100.)".format(stream, line),
		"NC_20_DELTAPHI_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_DELTAPHI', -100.)".format(stream, line),
		"NC_20_IT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_IT', -100.)".format(stream, line),
		"NC_20_MAXPT_PT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_MAXPT_PT', -100.)".format(stream, line),
		"NC_20_MAXPT_PX_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_MAXPT_PX', -100.)".format(stream, line),
		"NC_20_MAXPT_PY_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_MAXPT_PY', -100.)".format(stream, line),
		"NC_20_MAXPT_PZ_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'NC_MAXPT_PZ', -100.)".format(stream, line),
		"CCNC_20_IT_K" : "RELINFO('/Event/{0}/Phys/{1}/H_ConeIsoInfo_Cone20', 'CCNC_IT', -100.)".format(stream, line),

		# D+ (cone size = 2.0)
		"CC_20_ANGLE_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'CC_ANGLE', -100.)".format(stream, line),
		"CC_20_MULT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'CC_MULT', -100.)".format(stream, line),
		"CC_20_SPT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'CC_SPT', -100.)".format(stream, line),
		"CC_20_VPT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'CC_VPT', -100.)".format(stream, line),
		"CC_20_PX_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'CC_PX', -100.)".format(stream, line),
		"CC_20_PY_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'CC_PY', -100.)".format(stream, line),
		"CC_20_PZ_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'CC_PZ', -100.)".format(stream, line),
		"CC_20_PASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'CC_PASYM', -100.)".format(stream, line),
		"CC_20_PTASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'CC_PTASYM', -100.)".format(stream, line),
		"CC_20_PXASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'CC_PXASYM', -100.)".format(stream, line),
		"CC_20_PYASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'CC_PYASYM', -100.)".format(stream, line),
		"CC_20_PZASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'CC_PZASYM', -100.)".format(stream, line),
		"CC_20_DELTAETA_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'CC_DELTAETA', -100.)".format(stream, line),
		"CC_20_DELTAPHI_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'CC_DELTAPHI', -100.)".format(stream, line),
		"CC_20_PX_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'CC_PX', -100.)".format(stream, line),
		"CC_20_IT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'CC_IT', -100.)".format(stream, line),
		"CC_20_MAXPT_Q_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'CC_MAXPT_Q', -100.)".format(stream, line),
		"CC_20_MAXPT_PT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'CC_MAXPT_PT', -100.)".format(stream, line),
		"CC_20_MAXPT_PX_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'CC_MAXPT_PX', -100.)".format(stream, line),
		"CC_20_MAXPT_PY_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'CC_MAXPT_PY', -100.)".format(stream, line),
		"CC_20_MAXPT_PZ_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'CC_MAXPT_PZ', -100.)".format(stream, line),
		"CC_20_MAXPT_PE_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'CC_MAXPT_PE', -100.)".format(stream, line),
		"NC_20_ANGLE_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'NC_ANGLE', -100.)".format(stream, line),
		"NC_20_MULT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'NC_MULT', -100.)".format(stream, line),
		"NC_20_SPT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'NC_SPT', -100.)".format(stream, line),
		"NC_20_VPT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'NC_VPT', -100.)".format(stream, line),
		"NC_20_PX_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'NC_PX', -100.)".format(stream, line),
		"NC_20_PY_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'NC_PY', -100.)".format(stream, line),
		"NC_20_PZ_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'NC_PZ', -100.)".format(stream, line),
		"NC_20_PASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'NC_PASYM', -100.)".format(stream, line),
		"NC_20_PTASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'NC_PTASYM', -100.)".format(stream, line),
		"NC_20_PXASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'NC_PXASYM', -100.)".format(stream, line),
		"NC_20_PYASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'NC_PYASYM', -100.)".format(stream, line),
		"NC_20_PZASYM_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'NC_PZASYM', -100.)".format(stream, line),
		"NC_20_DELTAETA_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'NC_DELTAETA', -100.)".format(stream, line),
		"NC_20_DELTAPHI_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'NC_DELTAPHI', -100.)".format(stream, line),
		"NC_20_IT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'NC_IT', -100.)".format(stream, line),
		"NC_20_MAXPT_PT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'NC_MAXPT_PT', -100.)".format(stream, line),
		"NC_20_MAXPT_PX_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'NC_MAXPT_PX', -100.)".format(stream, line),
		"NC_20_MAXPT_PY_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'NC_MAXPT_PY', -100.)".format(stream, line),
		"NC_20_MAXPT_PZ_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'NC_MAXPT_PZ', -100.)".format(stream, line),
		"CCNC_20_IT_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_ConeIsoInfo_Cone20', 'CCNC_IT', -100.)".format(stream, line),

		# tau- (cone size = 2.0)
		"CC_20_ANGLE_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'CC_ANGLE', -100.)".format(stream, line),
		"CC_20_MULT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'CC_MULT', -100.)".format(stream, line),
		"CC_20_SPT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'CC_SPT', -100.)".format(stream, line),
		"CC_20_VPT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'CC_VPT', -100.)".format(stream, line),
		"CC_20_PX_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'CC_PX', -100.)".format(stream, line),
		"CC_20_PY_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'CC_PY', -100.)".format(stream, line),
		"CC_20_PZ_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'CC_PZ', -100.)".format(stream, line),
		"CC_20_PASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'CC_PASYM', -100.)".format(stream, line),
		"CC_20_PTASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'CC_PTASYM', -100.)".format(stream, line),
		"CC_20_PXASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'CC_PXASYM', -100.)".format(stream, line),
		"CC_20_PYASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'CC_PYASYM', -100.)".format(stream, line),
		"CC_20_PZASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'CC_PZASYM', -100.)".format(stream, line),
		"CC_20_DELTAETA_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'CC_DELTAETA', -100.)".format(stream, line),
		"CC_20_DELTAPHI_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'CC_DELTAPHI', -100.)".format(stream, line),
		"CC_20_PX_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'CC_PX', -100.)".format(stream, line),
		"CC_20_IT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'CC_IT', -100.)".format(stream, line),
		"CC_20_MAXPT_Q_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'CC_MAXPT_Q', -100.)".format(stream, line),
		"CC_20_MAXPT_PT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'CC_MAXPT_PT', -100.)".format(stream, line),
		"CC_20_MAXPT_PX_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'CC_MAXPT_PX', -100.)".format(stream, line),
		"CC_20_MAXPT_PY_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'CC_MAXPT_PY', -100.)".format(stream, line),
		"CC_20_MAXPT_PZ_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'CC_MAXPT_PZ', -100.)".format(stream, line),
		"CC_20_MAXPT_PE_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'CC_MAXPT_PE', -100.)".format(stream, line),
		"NC_20_ANGLE_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'NC_ANGLE', -100.)".format(stream, line),
		"NC_20_MULT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'NC_MULT', -100.)".format(stream, line),
		"NC_20_SPT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'NC_SPT', -100.)".format(stream, line),
		"NC_20_VPT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'NC_VPT', -100.)".format(stream, line),
		"NC_20_PX_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'NC_PX', -100.)".format(stream, line),
		"NC_20_PY_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'NC_PY', -100.)".format(stream, line),
		"NC_20_PZ_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'NC_PZ', -100.)".format(stream, line),
		"NC_20_PASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'NC_PASYM', -100.)".format(stream, line),
		"NC_20_PTASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'NC_PTASYM', -100.)".format(stream, line),
		"NC_20_PXASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'NC_PXASYM', -100.)".format(stream, line),
		"NC_20_PYASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'NC_PYASYM', -100.)".format(stream, line),
		"NC_20_PZASYM_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'NC_PZASYM', -100.)".format(stream, line),
		"NC_20_DELTAETA_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'NC_DELTAETA', -100.)".format(stream, line),
		"NC_20_DELTAPHI_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'NC_DELTAPHI', -100.)".format(stream, line),
		"NC_20_IT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'NC_IT', -100.)".format(stream, line),
		"NC_20_MAXPT_PT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'NC_MAXPT_PT', -100.)".format(stream, line),
		"NC_20_MAXPT_PX_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'NC_MAXPT_PX', -100.)".format(stream, line),
		"NC_20_MAXPT_PY_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'NC_MAXPT_PY', -100.)".format(stream, line),
		"NC_20_MAXPT_PZ_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'NC_MAXPT_PZ', -100.)".format(stream, line),
		"CCNC_20_IT_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_ConeIsoInfo_Cone20', 'CCNC_IT', -100.)".format(stream, line),

		# Vertex isolation BDT
		# B+
		"VTXISOBDTHARDFIRSTVALUE_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_VertexIsoBDTInfo', 'VTXISOBDTHARDFIRSTVALUE', -100.)".format(stream, line),
		"VTXISOBDTHARDSECONDVALUE_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_VertexIsoBDTInfo', 'VTXISOBDTHARDSECONDVALUE', -100.)".format(stream, line),
		"VTXISOBDTHARDTHIRDVALUE_B" : "RELINFO('/Event/{0}/Phys/{1}/BVars_VertexIsoBDTInfo', 'VTXISOBDTHARDTHIRDVALUE', -100.)".format(stream, line),

		# K+
		"VTXISOBDTHARDFIRSTVALUE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_VertexIsoBDTInfo', 'VTXISOBDTHARDFIRSTVALUE', -100.)".format(stream, line),
		"VTXISOBDTHARDSECONDVALUE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_VertexIsoBDTInfo', 'VTXISOBDTHARDSECONDVALUE', -100.)".format(stream, line),
		"VTXISOBDTHARDTHIRDVALUE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_VertexIsoBDTInfo', 'VTXISOBDTHARDTHIRDVALUE', -100.)".format(stream, line),

		# D+
		"VTXISOBDTHARDFIRSTVALUE_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_VertexIsoBDTInfo', 'VTXISOBDTHARDFIRSTVALUE', -100.)".format(stream, line),
		"VTXISOBDTHARDSECONDVALUE_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_VertexIsoBDTInfo', 'VTXISOBDTHARDSECONDVALUE', -100.)".format(stream, line),
		"VTXISOBDTHARDTHIRDVALUE_Dp" : "RELINFO('/Event/{0}/Phys/{1}/Dp_VertexIsoBDTInfo', 'VTXISOBDTHARDTHIRDVALUE', -100.)".format(stream, line),

		# D-
		"VTXISOBDTHARDFIRSTVALUE_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_VertexIsoBDTInfo', 'VTXISOBDTHARDFIRSTVALUE', -100.)".format(stream, line),
		"VTXISOBDTHARDSECONDVALUE_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_VertexIsoBDTInfo', 'VTXISOBDTHARDSECONDVALUE', -100.)".format(stream, line),
		"VTXISOBDTHARDTHIRDVALUE_Dm" : "RELINFO('/Event/{0}/Phys/{1}/Dm_VertexIsoBDTInfo', 'VTXISOBDTHARDTHIRDVALUE', -100.)".format(stream, line),

		# Track isolation BDT
		# K+
		"TRKISOBDTFIRSTVALUE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_TrackIsoBDTInfo', 'TRKISOBDTFIRSTVALUE', -100.)".format(stream, line),
		"TRKISOBDTSECONDVALUE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_TrackIsoBDTInfo', 'TRKISOBDTSECONDVALUE', -100.)".format(stream, line),
		"TRKISOBDTTHIRDVALUE_K" : "RELINFO('/Event/{0}/Phys/{1}/H_TrackIsoBDTInfo', 'TRKISOBDTTHIRDVALUE', -100.)".format(stream, line),

		# Dp K (D+ -> K- pi+ pi+)
		"TRKISOBDTFIRSTVALUE_Dp_K" : "RELINFO('/Event/{0}/Phys/{1}/Dp_pi2_TrackIsoBDTInfo', 'TRKISOBDTFIRSTVALUE', -100.)".format(stream, line),
		"TRKISOBDTSECONDVALUE_Dp_K" : "RELINFO('/Event/{0}/Phys/{1}/Dp_pi2_TrackIsoBDTInfo', 'TRKISOBDTSECONDVALUE', -100.)".format(stream, line),
		"TRKISOBDTTHIRDVALUE_Dp_K" : "RELINFO('/Event/{0}/Phys/{1}/Dp_pi2_TrackIsoBDTInfo', 'TRKISOBDTTHIRDVALUE', -100.)".format(stream, line),

		# Dp pi1
		"TRKISOBDTFIRSTVALUE_Dp_pi1" : "RELINFO('/Event/{0}/Phys/{1}/Dp_pi1_TrackIsoBDTInfo', 'TRKISOBDTFIRSTVALUE', -100.)".format(stream, line),
		"TRKISOBDTSECONDVALUE_Dp_pi1" : "RELINFO('/Event/{0}/Phys/{1}/Dp_pi1_TrackIsoBDTInfo', 'TRKISOBDTSECONDVALUE', -100.)".format(stream, line),
		"TRKISOBDTTHIRDVALUE_Dp_pi1" : "RELINFO('/Event/{0}/Phys/{1}/Dp_pi1_TrackIsoBDTInfo', 'TRKISOBDTTHIRDVALUE', -100.)".format(stream, line),

		# Dp pi2
		"TRKISOBDTFIRSTVALUE_Dp_pi2" : "RELINFO('/Event/{0}/Phys/{1}/Dp_pi3_TrackIsoBDTInfo', 'TRKISOBDTFIRSTVALUE', -100.)".format(stream, line),
		"TRKISOBDTSECONDVALUE_Dp_pi2" : "RELINFO('/Event/{0}/Phys/{1}/Dp_pi3_TrackIsoBDTInfo', 'TRKISOBDTSECONDVALUE', -100.)".format(stream, line),
		"TRKISOBDTTHIRDVALUE_Dp_pi2" : "RELINFO('/Event/{0}/Phys/{1}/Dp_pi3_TrackIsoBDTInfo', 'TRKISOBDTTHIRDVALUE', -100.)".format(stream, line),

		# Dm K (D- -> K+ pi- pi-)
		"TRKISOBDTFIRSTVALUE_Dm_K" : "RELINFO('/Event/{0}/Phys/{1}/Dm_pi2_TrackIsoBDTInfo', 'TRKISOBDTFIRSTVALUE', -100.)".format(stream, line),
		"TRKISOBDTSECONDVALUE_Dm_K" : "RELINFO('/Event/{0}/Phys/{1}/Dm_pi2_TrackIsoBDTInfo', 'TRKISOBDTSECONDVALUE', -100.)".format(stream, line),
		"TRKISOBDTTHIRDVALUE_Dm_K" : "RELINFO('/Event/{0}/Phys/{1}/Dm_pi2_TrackIsoBDTInfo', 'TRKISOBDTTHIRDVALUE', -100.)".format(stream, line),

		# Dm pi1
		"TRKISOBDTFIRSTVALUE_Dm_pi1" : "RELINFO('/Event/{0}/Phys/{1}/Dm_pi1_TrackIsoBDTInfo', 'TRKISOBDTFIRSTVALUE', -100.)".format(stream, line),
		"TRKISOBDTSECONDVALUE_Dm_pi1" : "RELINFO('/Event/{0}/Phys/{1}/Dm_pi1_TrackIsoBDTInfo', 'TRKISOBDTSECONDVALUE', -100.)".format(stream, line),
		"TRKISOBDTTHIRDVALUE_Dm_pi1" : "RELINFO('/Event/{0}/Phys/{1}/Dm_pi1_TrackIsoBDTInfo', 'TRKISOBDTTHIRDVALUE', -100.)".format(stream, line),

		# Dm pi2
		"TRKISOBDTFIRSTVALUE_Dm_pi2" : "RELINFO('/Event/{0}/Phys/{1}/Dm_pi3_TrackIsoBDTInfo', 'TRKISOBDTFIRSTVALUE', -100.)".format(stream, line),
		"TRKISOBDTSECONDVALUE_Dm_pi2" : "RELINFO('/Event/{0}/Phys/{1}/Dm_pi3_TrackIsoBDTInfo', 'TRKISOBDTSECONDVALUE', -100.)".format(stream, line),
		"TRKISOBDTTHIRDVALUE_Dm_pi2" : "RELINFO('/Event/{0}/Phys/{1}/Dm_pi3_TrackIsoBDTInfo', 'TRKISOBDTTHIRDVALUE', -100.)".format(stream, line),
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
	df = branch.addTupleTool("TupleToolDecayFit_DDK/df")
	df.constrainToOriginVertex = True

    # standard DTF
	dtf_DDK = branch.addTupleTool("TupleToolDecayTreeFitter/dtf")
	dtf_DDK.daughtersToConstrain = ['D+']
	dtf_DDK.UseFullTreeInName = True
	dtf_DDK.UpdateDaughters = True
	dtf_DDK.constrainToOriginVertex = True


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
    if isMC:
        addTupleToolL0Calo(DecayTreeTuple.Dp_K)
        addTupleToolL0Calo(DecayTreeTuple.Dp_pi1)
        addTupleToolL0Calo(DecayTreeTuple.Dp_pi2)
        addTupleToolL0Calo(DecayTreeTuple.Dm_K)
        addTupleToolL0Calo(DecayTreeTuple.Dm_pi1)
        addTupleToolL0Calo(DecayTreeTuple.Dm_pi2)
        addTupleToolL0Calo(DecayTreeTuple.Kp)

# Configure DaVinci
if isMC:
	DaVinci().TupleFile = 'DVntuple_MC_'+year+'_'+pol+'.root'
	# DaVinci().TupleFile = '/panfs/felician/B2DDK/ROOT_Sim/'+year+'/DVntuple_MC_'+year+'_'+pol+'.root'
else:
	DaVinci().TupleFile = 'DVntuple_data_'+year+'_'+pol+'.root'
	#DaVinci().TupleFile = '/panfs/felician/B2DDK/ROOT_Data/'+year+'/DVntuple_data_'+year+'_'+pol+'.root'

# Stream and stripping line we want to use
if isMC:
	stream = 'B2DDK.Strip' # for filtered MC
else:	
	stream = 'Bhadron'
line = 'B2KTauTau_DDKLine'
line_SS = 'B2KTauTau_DDKSSLine'

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

# Loose preselections to reduce the size of the tuples (in data)
# These (rectangular) selections are applied to MC offline
# filter code (will be updated)
filtercode = "( 2 == NINGENERATION( ('D+'==ABSID) & ( ADMASS('D+')<80*MeV ), 1 ) ) & ( 2 == NINGENERATION( ('K+'==ABSID) & ( PROBNNK > 0.1 ), 2 ) )"

if applySel:
	if isMC:
		b_strip = AutomaticData('/Event/{0}/Phys/{1}/Particles'.format(stream, line))
	if not isMC:
		b_strip    = AutomaticData('/Phys/{0}/Particles'.format(line)) 
		b_strip_SS = AutomaticData('/Phys/{0}/Particles'.format(line_SS)) 

	b_Filter    = FilterSelection("b_Filter", b_strip, Code = filtercode)
	dtt.Inputs = [b_Filter.outputLocation()]

	if not isMC:
		b_Filter_SS = FilterSelection("b_Filter_SS", b_strip_SS, Code = filtercode)
		dtt_SS.Inputs = [b_Filter_SS.outputLocation()]

# Do not apply rectangular cuts to MC @ DaVinci stage:
if not applySel:
	if isMC:
		dtt.Inputs = ['/Event/{0}/Phys/{1}/Particles'.format(stream, line)] 
	if not isMC:
		dtt.Inputs = ['Phys/{0}/Particles'.format(line)]
		dtt_SS.Inputs = ['Phys/{0}/Particles'.format(line_SS)]

dtt.setDescriptorTemplate('${Bp}[B+ -> ${Kp}K+ ${Dp}(D+ -> ${Dp_K}K- ${Dp_pi1}pi+ ${Dp_pi2}pi+) ${Dm}(D- -> ${Dm_K}K+ ${Dm_pi1}pi- ${Dm_pi2}pi-)]CC')
dtt_SS.setDescriptorTemplate('(${Bp}[B+ -> ${Kp}K+ ${Dp}(D+ -> ${Dp_K}K- ${Dp_pi1}pi+ ${Dp_pi2}pi+) ${Dm}(D+ -> ${Dm_K}K- ${Dm_pi1}pi+ ${Dm_pi2}pi+)]CC)'
						'||   (${Bp}[B+ -> ${Kp}K- ${Dp}(D+ -> ${Dp_K}K- ${Dp_pi1}pi+ ${Dp_pi2}pi+) ${Dm}(D+ -> ${Dm_K}K- ${Dm_pi1}pi+ ${Dm_pi2}pi+)]CC)'   )
# Generated (Acceptance cut: DaughtersInLHCb)
mcdtt.Decay = '[B+ => ^(D+ => ^K- ^pi+ ^pi+)  ^(D- => ^K+ ^pi- ^pi-) K+]CC' 

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

# Stripping filter
from PhysConf.Filters import LoKi_Filters
fltrs = LoKi_Filters (
    STRIP_Code = "HLT_PASS_RE('StrippingB2KTauTau_DDK(Line|SSLine)Decision')"
)
if not isMC: # don't apply to MC because it will affect the MCDecayTree
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
if not isMC:
	DaVinci().UserAlgorithms += [dtt_SS]

# Use the local input data
# from GaudiConf import IOHelper
# IOHelper('ROOT').inputFiles([
# 	'/panfs/felician/B2DDK/DST/00213826_00000018_1.b2ddk.strip.dst'
# ], clear=True)