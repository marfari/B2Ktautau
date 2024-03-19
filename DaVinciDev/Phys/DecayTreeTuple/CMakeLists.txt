###############################################################################
# (c) Copyright 2000-2019 CERN for the benefit of the LHCb Collaboration      #
#                                                                             #
# This software is distributed under the terms of the GNU General Public      #
# Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   #
#                                                                             #
# In applying this licence, CERN does not waive the privileges and immunities #
# granted to it by virtue of its status as an Intergovernmental Organization  #
# or submit itself to any jurisdiction.                                       #
###############################################################################
################################################################################
# Package: DecayTreeTuple
################################################################################
gaudi_subdir(DecayTreeTuple)

gaudi_depends_on_subdirs(Phys/DecayTreeFitter
                         Phys/MVADictTools
                         Phys/DecayTreeTupleBase
                         Phys/DecayTreeTupleDalitz
                         Phys/DecayTreeTupleJets
                         Phys/DecayTreeTupleMC
                         Phys/DecayTreeTupleMuonCalib
                         Phys/DecayTreeTupleReco
                         Phys/DecayTreeTupleTrigger
                         Phys/DecayTreeTupleTracking
                         Phys/DecayTreeTupleANNPID
                         Phys/TeslaTools
                         Phys/LoKiPhys)

find_package(HepMC)
find_package(Boost)
find_package(ROOT)
find_package(ROOT COMPONENTS Physics)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS})

gaudi_add_module(DecayTreeTuple
                 src/*.cpp
                 INCLUDE_DIRS HepMC
                 LINK_LIBRARIES HepMC DecayTreeFitter MVADictToolsLib DecayTreeTupleBaseLib LoKiPhysLib ROOT)

gaudi_install_python_modules()

gaudi_add_test(QMTest QMTEST)
