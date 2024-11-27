import pandas as pd
import xml.etree.ElementTree as ET
import numpy as np
import sys

def main(argv):

    year = argv[1]
    pol = argv[2]
    name = argv[3]

    if(name == "BuDDKp"):
        path = "/panfs/felician/SimulationJobs/12693500/201{0}/Sim10c_ReDecay/Mag{1}/xml/".format(year,pol)
    elif(name == "BdDDKp"):
        path = "/panfs/felician/SimulationJobs/11293500/201{0}/Sim10c_ReDecay/Mag{1}/xml/".format(year,pol)
    elif(name == "BsDDKp"):
        path = "/panfs/felician/SimulationJobs/13297500/201{0}/Sim10c_ReDecay/Mag{1}/xml/".format(year,pol)
    elif(name == "BuDDK0"):
        path = "/panfs/felician/SimulationJobs/12293500/201{0}/Sim10c_ReDecay/Mag{1}/xml/".format(year,pol)
    elif(name == "BdDDK0"):
        path = "/panfs/felician/SimulationJobs/11294100/201{0}/Sim10c_ReDecay/Mag{1}/xml/".format(year,pol)
    elif(name == "BuDD"):
        path = "/panfs/felician/SimulationJobs/12696700/201{0}/Sim10c_ReDecay/Mag{1}/xml/".format(year,pol)
    elif(name == "BdDD"):
        path = "/panfs/felician/SimulationJobs/11697700/201{0}/Sim10c_ReDecay/Mag{1}/xml/".format(year,pol)
    elif(name == "BsDD"):
        path = "/panfs/felician/SimulationJobs/13699600/201{0}/Sim10c_ReDecay/Mag{1}/xml/".format(year,pol)

    f = open(path+'201{0}_Mag{1}.txt'.format(year,pol))

    xml_data = f.readlines(0)
    xml = open(xml_data[0].rstrip('\n'), 'r').read()
    root = ET.XML(xml) 

    dataframes = []
    for i in range(len(xml_data)):
        data = []
        cols = []
        for i, child in enumerate(root):
            data.append([subchild.text for subchild in child])
            if(child.attrib == {'name': 'generator level cut'}):
                cols.append('efficiency')
                if cols[i] in cols:
                    cols[i] += '{0}'.format(i)
            else: 
                cols.append(child.tag)

        df = pd.DataFrame(data).T  # Write in DF and transpose it
        df.columns = cols  # Update column names
        dataframes.append(df)

    GEN_EVTs = 0
    ACC_EVTs = 0
    value = 0.
    error = 0.

    for i in range(len(dataframes)):
        dataframe = dataframes[i]
        ACC_EVTs += int(dataframe['efficiency'].iloc[0][-2])
        GEN_EVTs += int(dataframe['efficiency'].iloc[1][-2])
        value += float(dataframe['efficiency'].iloc[2][-2])
        error += float(dataframe['efficiency'].iloc[3][-2])

    print("Cocktail MC : ", name)
    print('ACC = ', ACC_EVTs, 'GEN = ', GEN_EVTs)
    #print('value = ', value/len(dataframes), 'error = ', error/len(dataframes))

if __name__ == "__main__":
    main(sys.argv)