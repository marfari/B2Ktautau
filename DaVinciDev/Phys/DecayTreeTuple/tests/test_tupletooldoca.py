###############################################################################
# (c) Copyright 2019 CERN for the benefit of the LHCb Collaboration           #
#                                                                             #
# This software is distributed under the terms of the GNU General Public      #
# Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   #
#                                                                             #
# In applying this licence, CERN does not waive the privileges and immunities #
# granted to it by virtue of its status as an Intergovernmental Organization  #
# or submit itself to any jurisdiction.                                       #
###############################################################################
from Configurables import (AnalysisConf, ApplicationMgr, DecayTreeTuple,
                           DstConf, LHCbApp, NTupleSvc, PhysConf, TrackSys)
from CommonParticles import StandardDstar
from DecayTreeTuple import Configuration
from PRConfig import TestFileDB

# Set up the DecayTreeTuple for D*+ -> (D0 -> X X) pi+
dtt = DecayTreeTuple("TupleDstToD0pi_D0ToHH")
dtt.Decay = "[D*(2010)+ -> X X]CC"
dtt.Inputs = list(StandardDstar.locations.keys())
dtt.ToolList = []
dtt.addBranches({"Dst": "[D*(2010)+ -> X X]CC"})

# Set up the tuple tool under test
docas = dtt.Dst.addTupleTool("TupleToolDOCA")
name, location1, location2 = list(
    zip(*[[
        "Dst_pi_D0_hm", "[D*(2010)+ -> (X -> X- X+) ^X+]CC",
        "[D*(2010)+ -> (X -> ^X- X+) X+]CC"
    ],
          [
              "Dst_pi_D0_hp", "[D*(2010)+ -> (X -> X- X+) ^X+]CC",
              "[D*(2010)+ -> (X -> X- ^X+) X+]CC"
          ],
          [
              "D0_hm_D0_hp", "[D*(2010)+ -> (X -> ^X- X+) X+]CC",
              "[D*(2010)+ -> (X -> X- ^X+) X+]CC"
          ]]))
docas.Name = list(name)
docas.Location1 = list(location1)
docas.Location2 = list(location2)

ApplicationMgr().TopAlg = [dtt]

app = LHCbApp()
app.EvtMax = 100

# Pick a file that has the reconstruction available
f = TestFileDB.test_file_db[
    "2017HLTCommissioning_BnoC_MC2015_Bs2PhiPhi_L0Processed"]
f.setqualifiers(configurable=app)
f.run(configurable=app, withDB=True)

print(app)

# Make sure underlying application configuration is run (in particular the
# DataOnDemand unpacking)
PhysConf().DataType = LHCbApp().DataType
AnalysisConf().DataType = LHCbApp().DataType
DstConf().DataType = LHCbApp().DataType
DstConf().EnableUnpack = PhysConf().EnableUnpack = [
    "Reconstruction", "Stripping"
]
TrackSys().DataType = LHCbApp().DataType

# Define the output ntuple name
ApplicationMgr().HistogramPersistency = "ROOT"
tuple_name = "test_tupletooldoca.root"
tuple_str = "FILE1 DATAFILE='{}' TYP='ROOT' OPT='NEW'".format(tuple_name)
ApplicationMgr().ExtSvc += [NTupleSvc(Output=[tuple_str])]
