##############################################################################
# (c) Copyright 2021 CERN for the benefit of the LHCb Collaboration           #
#                                                                             #
# This software is distributed under the terms of the GNU General Public      #
# Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   #
#                                                                             #
# In applying this licence, CERN does not waive the privileges and immunities #
# granted to it by virtue of its status as an Intergovernmental Organization  #
# or submit itself to any jurisdiction.                                       #
###############################################################################

import os
import sys

## START OF CONFIG

# Read comments and check vars at least until the end of config section

# List of input ROOT files with MC ntuples. Format:
#   (inputfile, outputfile, dataset)
# outputfile should be without the ".root" extension
# (it will be appended by the variable names when writing friend trees)
files = [
    ("root://myCerttxt@eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/m/mfaria/2024_07/909023/909023227/DVntuple_MC_Ktautau_2016_MagUp.root", "/panfs/felician/B2Ktautau/workflow/PIDCalib/Ktautau/2016_pidcorr", "MagUp_2016")
]

files_2016 = [
    ("root://myCerttxt@eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/m/mfaria/2024_07/909023/909023227/DVntuple_MC_Ktautau_2016_MagUp.root", "/panfs/felician/B2Ktautau/workflow/PIDCalib/Ktautau/2016_pidcorr", "MagUp_2016"),
    ("root://xrootd.echo.stfc.ac.uk//lhcb:user/lhcb/user/m/mfaria/2024_07/908789/908789714/DVntuple_MC_Ktautau_2016_MagUp.root", "/panfs/felician/B2Ktautau/workflow/PIDCalib/Ktautau/2016_pidcorr", "MagUp_2016"),
    ("root://myCerttxt@ccxrootdlhcb.in2p3.fr//pnfs/in2p3.fr/data/lhcb/LHCb_USER/lhcb/user/m/mfaria/2024_07/908789/908789744/DVntuple_MC_Ktautau_2016_MagUp.root", "/panfs/felician/B2Ktautau/workflow/PIDCalib/Ktautau/2016_pidcorr", "MagUp_2016"),
    ("root://myCerttxt@eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/m/mfaria/2024_07/908789/908789786/DVntuple_MC_Ktautau_2016_MagUp.root", "/panfs/felician/B2Ktautau/workflow/PIDCalib/Ktautau/2016_pidcorr", "MagUp_2016"),
    ("root://myCerttxt@ccxrootdlhcb.in2p3.fr//pnfs/in2p3.fr/data/lhcb/LHCb_USER/lhcb/user/m/mfaria/2024_07/908789/908789835/DVntuple_MC_Ktautau_2016_MagUp.root", "/panfs/felician/B2Ktautau/workflow/PIDCalib/Ktautau/2016_pidcorr", "MagUp_2016"),
    ("root://myCerttxt@eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/m/mfaria/2024_07/908789/908789809/DVntuple_MC_Ktautau_2016_MagUp.root", "/panfs/felician/B2Ktautau/workflow/PIDCalib/Ktautau/2016_pidcorr", "MagUp_2016"),
    ("root://myCerttxt@eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/m/mfaria/2024_07/908792/908792116/DVntuple_MC_Ktautau_2016_MagDown.root", "/panfs/felician/B2Ktautau/workflow/PIDCalib/Ktautau/2016_pidcorr", "MagDown_2016"),
    ("root://xrootd.echo.stfc.ac.uk//lhcb:user/lhcb/user/m/mfaria/2024_07/908792/908792126/DVntuple_MC_Ktautau_2016_MagDown.root", "/panfs/felician/B2Ktautau/workflow/PIDCalib/Ktautau/2016_pidcorr", "MagDown_2016"),
    ("root://myCerttxt@ccxrootdlhcb.in2p3.fr//pnfs/in2p3.fr/data/lhcb/LHCb_USER/lhcb/user/m/mfaria/2024_07/908792/908792139/DVntuple_MC_Ktautau_2016_MagDown.root", "/panfs/felician/B2Ktautau/workflow/PIDCalib/Ktautau/2016_pidcorr", "MagDown_2016"),
    ("root://myCerttxt@eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/m/mfaria/2024_07/908792/908792153/DVntuple_MC_Ktautau_2016_MagDown.root", "/panfs/felician/B2Ktautau/workflow/PIDCalib/Ktautau/2016_pidcorr", "MagDown_2016"),
    ("root://xrootd-lhcb.cr.cnaf.infn.it:1094//storage/gpfs_lhcb/lhcb/user/m/mfaria/2024_07/908792/908792168/DVntuple_MC_Ktautau_2016_MagDown.root", "/panfs/felician/B2Ktautau/workflow/PIDCalib/Ktautau/2016_pidcorr", "MagDown_2016")
]

files_2017 = [
    ("root://myCerttxt@eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/m/mfaria/2024_07/908794/908794696/DVntuple_MC_Ktautau_2017_MagUp.root", "/panfs/felician/B2Ktautau/workflow/PIDCalib/Ktautau/2017_pidcorr", "MagUp_2017"),
    ("root://xrootd.echo.stfc.ac.uk//lhcb:user/lhcb/user/m/mfaria/2024_07/908794/908794699/DVntuple_MC_Ktautau_2017_MagUp.root", "/panfs/felician/B2Ktautau/workflow/PIDCalib/Ktautau/2017_pidcorr", "MagUp_2017"),
    ("root://myCerttxt@ccxrootdlhcb.in2p3.fr//pnfs/in2p3.fr/data/lhcb/LHCb_USER/lhcb/user/m/mfaria/2024_07/908794/908794703/DVntuple_MC_Ktautau_2017_MagUp.root", "/panfs/felician/B2Ktautau/workflow/PIDCalib/Ktautau/2017_pidcorr", "MagUp_2017"),
    ("root://myCerttxt@eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/m/mfaria/2024_07/908795/908795152/DVntuple_MC_Ktautau_2017_MagDown.root", "/panfs/felician/B2Ktautau/workflow/PIDCalib/Ktautau/2017_pidcorr", "MagDown_2017"),
    ("root://myCerttxt@ccxrootdlhcb.in2p3.fr//pnfs/in2p3.fr/data/lhcb/LHCb_USER/lhcb/user/m/mfaria/2024_07/908795/908795170/DVntuple_MC_Ktautau_2017_MagDown.root", "/panfs/felician/B2Ktautau/workflow/PIDCalib/Ktautau/2017_pidcorr", "MagDown_2017"),
    ("root://xrootd.echo.stfc.ac.uk//lhcb:user/lhcb/user/m/mfaria/2024_07/908795/908795186/DVntuple_MC_Ktautau_2017_MagDown.root", "/panfs/felician/B2Ktautau/workflow/PIDCalib/Ktautau/2017_pidcorr", "MagDown_2017")
]

files_2018 = [
    ("root://myCerttxt@ccxrootdlhcb.in2p3.fr//pnfs/in2p3.fr/data/lhcb/LHCb_USER/lhcb/user/m/mfaria/2024_07/908796/908796765/DVntuple_MC_Ktautau_2018_MagUp.root", "/panfs/felician/B2Ktautau/workflow/PIDCalib/Ktautau/2018_pidcorr", "MagUp_2018"),
    ("root://myCerttxt@eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/m/mfaria/2024_07/908796/908796776/DVntuple_MC_Ktautau_2018_MagUp.root", "/panfs/felician/B2Ktautau/workflow/PIDCalib/Ktautau/2018_pidcorr", "MagUp_2018"),
    ("root://myCerttxt@eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/m/mfaria/2024_07/908796/908796784/DVntuple_MC_Ktautau_2018_MagUp.root", "/panfs/felician/B2Ktautau/workflow/PIDCalib/Ktautau/2018_pidcorr", "MagUp_2018"),
    ("root://myCerttxt@lhcbdcache-kit.gridka.de:1094//pnfs/gridka.de/lhcb/LHCb_USER/lhcb/user/m/mfaria/2024_07/908796/908796793/DVntuple_MC_Ktautau_2018_MagUp.root", "/panfs/felician/B2Ktautau/workflow/PIDCalib/Ktautau/2018_pidcorr", "MagUp_2018"),
    ("root://myCerttxt@lhcbdcache-kit.gridka.de:1094//pnfs/gridka.de/lhcb/LHCb_USER/lhcb/user/m/mfaria/2024_07/908796/908796808/DVntuple_MC_Ktautau_2018_MagUp.root", "/panfs/felician/B2Ktautau/workflow/PIDCalib/Ktautau/2018_pidcorr", "MagUp_2018"),
    ("root://xrootd.echo.stfc.ac.uk//lhcb:user/lhcb/user/m/mfaria/2024_07/908796/908796799/DVntuple_MC_Ktautau_2018_MagUp.root", "/panfs/felician/B2Ktautau/workflow/PIDCalib/Ktautau/2018_pidcorr", "MagUp_2018"),
    ("root://myCerttxt@ccxrootdlhcb.in2p3.fr//pnfs/in2p3.fr/data/lhcb/LHCb_USER/lhcb/user/m/mfaria/2024_07/908799/908799635/DVntuple_MC_Ktautau_2018_MagDown.root", "/panfs/felician/B2Ktautau/workflow/PIDCalib/Ktautau/2018_pidcorr", "MagDown_2018"),
    ("root://myCerttxt@eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/m/mfaria/2024_07/908799/908799651/DVntuple_MC_Ktautau_2018_MagDown.root", "/panfs/felician/B2Ktautau/workflow/PIDCalib/Ktautau/2018_pidcorr", "MagDown_2018"),
    ("root://myCerttxt@eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/m/mfaria/2024_07/909023/909023177/DVntuple_MC_Ktautau_2018_MagDown.root", "/panfs/felician/B2Ktautau/workflow/PIDCalib/Ktautau/2018_pidcorr", "MagDown_2018"),
    ("root://myCerttxt@eoslhcb.cern.ch//eos/lhcb/grid/user/lhcb/user/m/mfaria/2024_07/908799/908799672/DVntuple_MC_Ktautau_2018_MagDown.root", "/panfs/felician/B2Ktautau/workflow/PIDCalib/Ktautau/2018_pidcorr", "MagDown_2018"),
    ("root://xrootd.echo.stfc.ac.uk//lhcb:user/lhcb/user/m/mfaria/2024_07/908799/908799690/DVntuple_MC_Ktautau_2018_MagDown.root", "/panfs/felician/B2Ktautau/workflow/PIDCalib/Ktautau/2018_pidcorr", "MagDown_2018"),
    ("root://myCerttxt@xrootd.grid.surfsara.nl//pnfs/grid.sara.nl/data/lhcb/LHCb_USER/lhcb/user/m/mfaria/2024_07/908799/908799684/DVntuple_MC_Ktautau_2018_MagDown.root", "/panfs/felician/B2Ktautau/workflow/PIDCalib/Ktautau/2018_pidcorr", "MagDown_2018")
]

files = files_2016

# Name of the input tree
# Could also include ROOT directory, e.g. "Dir/Tree"
input_tree = "ntuple/DecayTree"

# Postfixes of the Pt, Eta and Ntracks variables (ntuple variable name the particle name part)
# e.g. if the ntuple contains "pi_PT" for the particle "pi", it should be just "PT"
ptvar  = "PT"
etavar = "ETA"
pvar   = None
## Could also use P variable instead of eta
#etavar = None
#pvar   = "p"
ntrvar = "nTracks"  # Number of tracks variable, this should correspond to the number of "Best tracks", not "Long tracks"!

seed = 10000        # initial random seed for resampling
friend = False      # If True, write friend trees with resampled PID instead of copying the whole tree

scale = (1, 1, 1.15)  # List of scale factors for each dimension of input data, to increase track multiplicity by 15%

# List of (kernel, seeds) combinations for raw template smearing
kernels = [
  ('default',    [0, 1, 2, 3]),   # Default kernel width, no template sampling 
                                  # and three "bootstrapped" templates with seeds 1..3
  ('syst1',      [0]),            # Wider "syst1" kernel, no template sampling 
  ((2, 2, 2, 2), [0])             # Narrower kernel with custom widths (2,2,2,2) bins in PID,Pt,Eta,Ntracks 
                                  # (will appear as "kern2" in ntuple), no template sampling
]

# Configuration dictionary for resampling, in the form {particle_name}:{pidvars}
# For each {particle_name}, {pidvars} is a dictionary in the form {ntuple_variable}:({sample}, {PID_var}, {kernels}),
#   where
#     {ntuple_variable} is the name of the corresponding ntuple PID variable without the particle name part
#                       (e.g. for "pi_PIDK" branch of particle "pi" it should be just "PIDK"); 
#     {sample} is the name of the calibration sample
#     {PID_var} is the string describing the PID variable template.
#     {kernels} is the list of kernels for template smearing (see example above)
# 
# Run pidgen2.list_variables to get the full list of PID variables
# Run pidgen2.list_samples to get the full list of calibration samples

config = {
  'Kp' :    {
            "ProbNNk"    : ("K_Dstar2Dpi", "MC15TuneV1_ProbNNK", kernels),
            "ProbNNpi"   : ("K_Dstar2Dpi", "MC15TuneV1_ProbNNpi", kernels)
           },
  'taup_pi1' :    {
            "ProbNNk"    : ("pi_Dstar2Dpi", "MC15TuneV1_ProbNNK", kernels),
            "ProbNNpi"   : ("pi_Dstar2Dpi", "MC15TuneV1_ProbNNpi", kernels)
           },
  'taup_pi2' :    {
            "ProbNNk"    : ("pi_Dstar2Dpi", "MC15TuneV1_ProbNNK", kernels),
            "ProbNNpi"   : ("pi_Dstar2Dpi", "MC15TuneV1_ProbNNpi", kernels)
           },
  'taup_pi3' :    {
            "ProbNNk"    : ("pi_Dstar2Dpi", "MC15TuneV1_ProbNNK", kernels),
            "ProbNNpi"   : ("pi_Dstar2Dpi", "MC15TuneV1_ProbNNpi", kernels)
           },
  'taum_pi1' :    {
            "ProbNNk"    : ("pi_Dstar2Dpi", "MC15TuneV1_ProbNNK", kernels),
            "ProbNNpi"   : ("pi_Dstar2Dpi", "MC15TuneV1_ProbNNpi", kernels)
           },
  'taum_pi2' :    {
            "ProbNNk"    : ("pi_Dstar2Dpi", "MC15TuneV1_ProbNNK", kernels),
            "ProbNNpi"   : ("pi_Dstar2Dpi", "MC15TuneV1_ProbNNpi", kernels)
           },
  'taum_pi3' :    {
            "ProbNNk"    : ("pi_Dstar2Dpi", "MC15TuneV1_ProbNNK", kernels),
            "ProbNNpi"   : ("pi_Dstar2Dpi", "MC15TuneV1_ProbNNpi", kernels)
           }
}

## END OF CONFIG

output_tree = input_tree.split("/")[-1]
treename = input_tree

from pidgen2.resample import resample

varseed = seed

for input_file, output_file, dataset in files :

  input_location = f"{input_file}:{input_tree}"

  for track, subst in config.items() :
    for var, (sample, calibvar, kernel_list) in subst.items() :

      # Create the list of input branches, depending on whether Eta or P variable is available
      if pvar is None : 
        branches = f"{track}_{ptvar}:{track}_{etavar}:{ntrvar}"
        eta_from_p = False
      else : 
        branches = f"{track}_{ptvar}:{track}_{pvar}:{ntrvar}"
        eta_from_p = True

      if friend : 
         output_root_file = f"{output_file}_{track}_{var}.root"
      else : 
         output_root_file = f"{output_file}.root"

      varseed += 1   # Increment random seed to make sure all resampled variables are independent

      # Run resampling of a single variable in a single file
      resample(
         input = input_location,    # Input tuple
         sample = sample,           # Calibration sample (e.g. "pi_Dstar2Dpi")
         dataset = dataset,         # Dataset (e.g. "MagUp_2016")
         variable = calibvar,       # Calibration variable (e.g. "MC15TuneV1_ProbNNK")
         branches = branches,       # List of resampling branches (typically, corresponding to Pt, Eta and Ntracks, e.g. "pt:eta:ntr")
         output = output_root_file, # Output ROOT file name
         outtree = output_tree,     # Output tree name
         plot = True,               # If template needs to be created from scratch, produce control plots
         pidgen = f"{track}_{var}_pidgen",   # Name of the resampled PID branch
         stat = f"{track}_{var}_pidstat",    # Name of output branch with calibration statistics for each resampled event
         resampling_seed = varseed, # Random seed for resampling
         kernels = kernel_list,     # List of kernels and template seeds
         verbose = False,           # Print debugging information
         eta_from_p = eta_from_p,   # If eta needs to be calculated from p and pt
         friend = friend,           # If True, write output to friend trees
         library = "ak",            # Library to handle ROOT files with uproot if friend=False, can be 
                                    # "ak" (Awkward array), "np" (numpy) or "pd" (Pandas)
         step_size = 50000,         # Chunk size when writing files with uproot if friend=False
         nan = -1000.,              # Numerical value to substitute NaN, for regions w/o calibration data 
         scale = scale,             # Scale factors for input data 
         local_storage = "./templates/",  # Directory for local template storage, used when the template is not available in 
                                          # the global storage on EOS. 
      )

      if not friend : 
        # All subsequent calls use output file as input
        input_location = f"{output_root_file}:{output_tree}"