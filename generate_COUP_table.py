import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from astropy import units as u
import matplotlib as mlp
from astropy.constants import h, k_B, c

distancePcTaurus = 140. * u.pc
distanceCMTaurus = distancePcTaurus.to(u.cm).value

distancePcOrion = 450. * u.pc
distanceCMOrion = distancePcOrion.to(u.cm).value

def generate_COUP_table(modelfilename,lumfilename):
    '''Takes data from two separate tables in Getman et al. (2005) and selects needed columns from them both for use in plotting HL Tau 
    and XZ Tau data against them.'''
    #Read tables in
    tableCOUPModels = Table.read(modelfilename,format='mrt')
    tableCOUPLuminosities = Table.read(lumfilename,format='mrt')

    #Generate new columns for absorbed and unabsorbed luminosity, and absorbed and unabsorbed flux.
    tableCOUPLuminosities['LumAbs'] = 0.
    tableCOUPLuminosities['LumUnabs'] = 0.
    tableCOUPLuminosities['FluxAbs'] = 0.
    tableCOUPLuminosities['FluxUnabs'] = 0.

    #Populate new columns from data in tables
    for i in range(len(tableCOUPLuminosities['Seq'])):
        tableCOUPLuminosities['LumAbs'][i]     = 10.**(tableCOUPLuminosities['Lt'][i])
        tableCOUPLuminosities['LumUnabs'][i]   = 10.**(tableCOUPLuminosities['Ltc'][i])
        tableCOUPLuminosities['FluxAbs'][i]    = tableCOUPLuminosities['LumAbs'][i] / (4. * np.pi * (distanceCMOrion**2))
        tableCOUPLuminosities['FluxUnabs'][i]  = tableCOUPLuminosities['LumAbs'][i] / (4. * np.pi * (distanceCMOrion**2))

    #Map the rows based on label for each table to ensure that the correct rows go with the correct object from each table.
    mapIndexLuminosities = {tableCOUPLuminosities['CXOONCJ'][i]: i for i in range(len(tableCOUPLuminosities['Seq']))}
    mapIndexModels       = {tableCOUPModels['CXOONCJ'][i]: i for i in range(len(tableCOUPModels['Seq']))}

    #Add new columns to the Models table to populate with info from the Luminosities table
    tableCOUPModels['LogLum'] = 0.
    tableCOUPModels['FluxAbs'] = 0.
    tableCOUPModels['FluxUnabs'] = 0.
    tableCOUPModels['LogFluxUnabs'] = 0.
    tableCOUPModels['EMRatio'] = 0.
    tableCOUPModels['LogEMRatio'] = 0.
    tableCOUPModels['LogkT1'] = 0.
    tableCOUPModels['LogkT2'] = 0.

    #Populate the new columns
    for i in range(len(tableCOUPModels['Seq'])):
        CXkey = tableCOUPModels['CXOONCJ'][i]
        lumindex = mapIndexLuminosities[CXkey]
        
        tableCOUPModels['LogLum'][i] = tableCOUPLuminosities['Ltc'][lumindex]
        tableCOUPModels['FluxAbs'][i] = tableCOUPLuminosities['FluxAbs'][lumindex]
        tableCOUPModels['FluxUnabs'][i] = tableCOUPLuminosities['FluxUnabs'][lumindex]
        tableCOUPModels['LogFluxUnabs'][i] = np.log10(tableCOUPModels['FluxUnabs'][lumindex])
        tableCOUPModels['EMRatio'][i] = np.log10(tableCOUPModels['EMRatio'][i])
        tableCOUPModels['LogkT1'][i]  = np.log10(tableCOUPModels['kT1'][i])
        tableCOUPModels['LogkT2'][i]  = np.log10(tableCOUPModels['kT2'][i])

    #Identify which objects in each table to use--no upper limits in luminosity, no 'm' or 'p' flags in models, and valid comparisons.
    CXkeysUseLum   = [tableCOUPLuminosities['CXOONCJ'][i] for i in range(len(tableCOUPLuminosities['CXOONCJ'])) 
                      if ((not ('<' in tableCOUPLuminosities['l_Ls'][i])) and (not ('<' in tableCOUPLuminosities['l_Lh'][i])))]
    
    CXkeysUseMod   = [tableCOUPModels['CXOONCJ'][i] for i in range(len(tableCOUPModels['CXOONCJ'])) 
                      if ((not ('m' in tableCOUPModels['FitFl'][i])) and (not ('p' in tableCOUPModels['FitFl'][i])))]
    
    CXkeysUseComps = [tableCOUPModels['CXOONCJ'][i] for i in range(len(tableCOUPModels['CXOONCJ'])) 
                      if ((not (type(tableCOUPModels['kT2'][i])) == (type(tableCOUPModels['kT2'][2]))) and 
                          (tableCOUPModels['kT2'][i] < 15.))]
    
    CXkeysUse = [x for x in CXkeysUseLum if (x in CXkeysUseMod) and (x in CXkeysUseComps)]

    #return the table with good objects
    return tableCOUPModels[[mapIndexModels[x] for x in CXkeysUse]]