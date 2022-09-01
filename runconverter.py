from Gaudi.Configuration import *

from Configurables import DaVinci

from Configurables import FilterDesktop

from StandardParticles import StdAllNoPIDsPions

#from Configurables import JsonConverter

myconverter = JsonConverter()

DaVinci().EventPreFilters = [ ]

DaVinci().UserAlgorithms = [ myconverter ]

DaVinci().Input = ['/panfs/felician/B2Ktautau/DST_Data/MagDown/2018/00092561_00000274_1.bhadron.mdst']

DaVinci().DataType = "2018"
DaVinci().Simulation = False

DaVinci().PrintFreq = 1

#DaVinci().DDDBtag = "dddb-20130111"
#DaVinci().CondDBtag = "cond-20130111"
DaVinci().CondDBtag = "cond-20140604"
DaVinci().DDDBtag = "dddb-20130929-1"