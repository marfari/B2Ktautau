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
gaudi_subdir(DecayTreeFitter)

gaudi_depends_on_subdirs(Phys/DaVinciKernel
                         Phys/LoKiCore
                         Tr/TrackInterfaces
                         Tr/TrackKernel)

find_package(Boost)
find_package(ROOT)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS})

# when using Gaudi with new CMake configuration CLHEP has to be used in
# a slightly different way
if (NOT TARGET CLHEP::CLHEP)
  find_package(CLHEP)
  include_directories(SYSTEM ${CLHEP_INCLUDE_DIRS})
  set(clhep_inc CLHEP)
  set(clhep_link CLHEP)
else()
  set(clhep_inc)
  set(clhep_link CLHEP::CLHEP)
endif()

gaudi_add_library(DecayTreeFitter
                  src/*.cpp
                  PUBLIC_HEADERS DecayTreeFitter
                  INCLUDE_DIRS ${clhep_inc} AIDA Tr/TrackInterfaces
                  LINK_LIBRARIES ${clhep_link} DaVinciKernelLib LoKiCoreLib TrackKernel)

gaudi_add_dictionary(DecayTreeFitter
                     dict/DecayTreeFitterDict.h
                     dict/DecayTreeFitterDict.xml
                     INCLUDE_DIRS Tr/TrackInterfaces
                     LINK_LIBRARIES DaVinciKernelLib LoKiCoreLib TrackKernel DecayTreeFitter)

