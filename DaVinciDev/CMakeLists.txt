CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)

#---------------------------------------------------------------
# Load macros and functions for Gaudi-based projects
find_package(GaudiProject)
#---------------------------------------------------------------

if(EXISTS ${CMAKE_SOURCE_DIR}/Online/RootCnv)
  set(override_subdirs RootCnv)
endif()

# Declare project name and version
# Syntax:
#    gaudi_project(this_project this_version
#                  USE dep_project_1 version_1 [project_2 version_2 ...]
#                  [DATA pkg1 [VERSION vers1] [pkg2 ...]])
gaudi_project(DaVinciDev HEAD
              USE DaVinci v45r6)

# Enable Ganga integration
include(/cvmfs/lhcb.cern.ch/lib/var/lib/LbEnv/2399/stable/linux-64/lib/python3.9/site-packages/LbDevTools/data/cmake/GangaTools.cmake)
enable_ganga_integration()
