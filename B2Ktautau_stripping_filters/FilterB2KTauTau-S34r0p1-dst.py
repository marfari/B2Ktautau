"""
Stripping Filtering file for B2KTauTau S34r0p1 line, DST output (2018)
@author Maria Carolina Feliciano Faria
@date   2022-12-15
"""

#stripping version
stripping='stripping34r0p1'

#use CommonParticlesArchive
from CommonParticlesArchive import CommonParticlesArchiveConf
CommonParticlesArchiveConf().redirect(stripping)

from Gaudi.Configuration import *
MessageSvc().Format = "% F%30W%S%7W%R%T %0W%M"

# Disable the cache in Tr/TrackExtrapolators
from Configurables import TrackStateProvider
TrackStateProvider().CacheStatesOnDemand = False

#Fix for TrackEff lines
from Configurables import DecodeRawEvent
DecodeRawEvent().setProp("OverrideInputs",4.2)

# Build the streams and stripping object
from StrippingConf.Configuration import StrippingConf, StrippingStream
from StrippingSettings.Utils import strippingConfiguration
from StrippingArchive.Utils import buildStreams
from StrippingArchive import strippingArchive

#get the configuration dictionary from the database
config  = strippingConfiguration(stripping)
#get the line builders from the archive
archive = strippingArchive(stripping)
streams = buildStreams(stripping = config, archive = archive)

# Merge into one stream and run in filter mode
AllStreams = StrippingStream("B2KTauTau.Strip")

linesToAdd = []
for stream in streams:
    for line in stream.lines:
        if 'B2KTauTauLine' in line.name():
            line._prescale = 1.0
            linesToAdd.append(line)
AllStreams.appendLines(linesToAdd)

sc = StrippingConf( Streams = [ AllStreams ],
                    MaxCandidates = 2000,
                    TESPrefix = 'Strip'
                    )

AllStreams.sequence().IgnoreFilterPassed = False # so we filter events by stripping selected lines

# Configure the dst writers for the output
enablePacking = True
#enablePacking = False

from DSTWriters.microdstelements import *
from DSTWriters.Configuration import ( SelDSTWriter,
                                       stripDSTStreamConf,
                                       stripDSTElements
                                       )

#############################################################################
# Configuration of SelDSTWriter
SelDSTWriterElements = {
    'default'              : stripDSTElements(stripPrefix='Strip')
    }


SelDSTWriterConf = {
    'default'              : stripDSTStreamConf(stripPrefix='Strip')
    }

for stream in sc.activeStreams() :
   print "there is a stream called " + stream.name() + " active"

dstWriter = SelDSTWriter( "MyDSTWriter",
                          StreamConf = SelDSTWriterConf,
                          MicroDSTElements = SelDSTWriterElements,
                          OutputFileSuffix ='RD',
                          SelectionSequences = sc.activeStreams()
                          )

# Add stripping TCK 
from Configurables import StrippingTCK
stck = StrippingTCK(HDRLocation = '/Event/Strip/Phys/DecReports', TCK=0x44103401)
#0xVVVVSSSS, VVVV DaVinci version, SSSS Stripping number

# DaVinci Configuration
from Configurables import DaVinci
DaVinci().EvtMax = -1 # Number of events
DaVinci().Simulation = True
DaVinci().HistogramFile = "DVHistos.root"
DaVinci().appendToMainSequence( [ sc.sequence() ] )
DaVinci().appendToMainSequence( [ stck ] )
DaVinci().appendToMainSequence( [ dstWriter.sequence() ] )
DaVinci().ProductionType = "Stripping"
#DaVinci().DataType = "2016"

# Change the column size of Timing table
from Configurables import TimingAuditor, SequencerTimerTool
TimingAuditor().addTool(SequencerTimerTool,name="TIMER")
TimingAuditor().TIMER.NameSize = 60

# Save info of rejection factor in the filtering
#from Configurables import DumpFSR
#DumpFSR().OutputLevel = 3
#DumpFSR().AsciiFileName = "dumpfsr.check.output.txt"
#DaVinci().MoniSequence += [ DumpFSR() ]
