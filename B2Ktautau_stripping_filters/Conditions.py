from Configurables import DaVinci

isMC = True
year = '2016'
pol = 'MagDown'

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