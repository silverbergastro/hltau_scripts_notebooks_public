import numpy as np
from ciao_contrib.runtool import *
from sherpa.astro.ui import *
import matplotlib.pyplot as plt
import copy
from astropy.table import Table
from astropy import units as u
import time
import pprint

def set_model_from_dict(inputdict):
    '''Set the parameters for a Sherpa plasma model defined with components 'a1', 's1', and 's2.' Uses the 'setattr' capability to set the value based on the provided dictionary, using keysplit to divide the dictionary key into model component and model attribute.'''
    for key in inputdict.keys():
        if (key[-1] in ['+-']) or ('.' not in key):
            continue
        keysplit = key.split('.')
        if keysplit[0] == 'a1':
            setattr(a1,keysplit[1],inputdict[key])
        if keysplit[0] == 's1':
            setattr(s1,keysplit[1],inputdict[key])
        if keysplit[0] == 's2':
            setattr(s2,keysplit[1],inputdict[key])
    return


def set_new_groupings(obsid,binsize,lowlim,hilim):
    '''A simple function for re-defining the grouping of counts for a given observation.'''
    ungroup(obsid)
    notice(obsid)
    notice_id(obsid,lowlim,hilim)
    group_counts(obsid,binsize,tabStops=~(get_data(obsid).mask))
    return


def get_results_dict(obsidlabel,results):
    '''Write best-fit parameters from Sherpa's 'fit' command to a dictionary for future retrieval.'''
    resultdict = {'obsid':obsidlabel,'rstat': results.rstat}
    for i in range(len(results.parnames)):
        resultdict[results.parnames[i]] = results.parvals[i] 
    return resultdict

    
def get_confresults_dict(resultsdict, confresults):
    '''Write results from Sherpa's 'conf' command (best-fit parameters and uncertainties around those parameters) to a dictionary.'''
    confresultsdict = {x: resultsdict[x] for x in resultsdict.keys() if '.' not in x}
    for i in range(len(confresults.parnames)):
        confresultsdict[confresults.parnames[i]]     = confresults.parvals[i]
        confresultsdict[confresults.parnames[i]+'-'] = confresults.parmins[i]
        confresultsdict[confresults.parnames[i]+'+'] = confresults.parmaxes[i]
    return confresultsdict



def plot_chandra_data(baseid,baseidUse,filestem,gratings=False,binning=0,useinputbins=False):
    '''Generate model and data plots of HL Tau and write to file. If gratings are used, generate plots for zeroth order and gratings.'''

    #Plot in logspace, with energy as the x-axis
    set_xlog()
    set_ylog()
    set_analysis(baseidUse,'energy')

    print(baseid,baseidUse)
    if not useinputbins:
        #Re-bin to 15 counts for cleaner plotting
        for specid in [baseidUse]:
            ungroup(specid)
            notice_id(specid,0.5,8.)
            mask = get_data(specid).mask
            group_counts(specid,15,tabStops=~mask)
    plot_fit_delchi(baseidUse,color='C0')
    plt.tight_layout()
    plt.title(str(baseid))
    plt.savefig('/data/swolk/SILVERBERG/hltau_plots/'+str(baseid)+filestem)
    #Reset binning for future analysis
    if not useinputbins:
        for specid in [baseidUse]:
            ungroup(specid)
            if binning > 0:
                notice_id(specid,0.5,8.)
                mask = get_data(specid).mask
                group_counts(specid,binning,tabStops=~mask)
    
    if gratings:
        #Plot gratings in linear space with wavelength as the x-axis
        set_xlinear()
        set_ylinear()
        set_analysis(baseidUse+1,'wavelength')
        set_analysis(baseidUse+2,'wavelength')
        print(baseid,baseidUse+1,baseidUse+2)
        if not useinputbins:
            for specid in [baseidUse+1, baseidUse+2]:
                ungroup(specid)
                mask = get_data(specid).mask
                group_adapt(specid,2,tabStops=~mask)
        
        plot_fit_delchi(baseidUse+1,color='C1')
        plot_fit_delchi(baseidUse+2,color='C2',overplot=True)
        
        [ax1,ax2] = plt.gcf().axes
        print(ax1.lines)
        ax1.lines[0].set_linestyle('-')
        ax1.lines[1].set_marker('')
        ax1.lines[2].set_linestyle('-')
        ax1.lines[3].set_marker('')
        ax1.lines[4].set_linestyle('-')
        ax1.lines[1].set_marker('')
        ax1.lines[2].set_linestyle('-')
        ax1.lines[3].set_marker('')
        ax1.set_xlim([1.,12.])
        ax2.set_xlim([1.,12.])
        plt.tight_layout()
        plt.title(str(baseid))
        plt.savefig('/data/swolk/SILVERBERG/hltau_plots/'+str(baseid)+'_gratings'+filestem)
        #Reset for future analysis
        set_analysis(baseidUse+1,'energy')
        set_analysis(baseidUse+2,'energy')
        if not useinputbins:
            for specid in [baseidUse+1, baseidUse+2]:
                ungroup(specid)
                if binning > 0:
                    notice_id(specid,0.5,8.)
                    mask = get_data(specid).mask
                    group_counts(specid,binning,tabStops=~mask)
    return

def plot_xmm_data(baseid,baseidUse,filestem,plotMOS2=True,binning=0,useinputbins=False):
    '''Plot XMM data, adapting plotting binning for more readable plots and taking into account which detectors are used in each plot.'''
    set_xlog()
    set_ylog()
    set_analysis(baseidUse+1,'energy')
    set_analysis(baseidUse+2,'energy')
    set_analysis(baseidUse+3,'energy')
    print(baseid,baseidUse)

    #Re-bin for plotting
    if not useinputbins:
        for specid in [baseidUse+1,baseidUse+2,baseidUse+3]:
            ungroup(specid)
            notice_id(specid,0.3,9.)
            mask = get_data(specid).mask
            group_counts(specid,15,tabStops=~mask)
    plot_fit_delchi(baseidUse+1,color='C0')
    #Only plot MOS2 if MOS2 was used in the fit
    if plotMOS2:
        plot_fit_delchi(baseidUse+2,color='C1',overplot=True)
    plot_fit_delchi(baseidUse+3,color='C2',overplot=True)
    plt.tight_layout()
    plt.title(str(baseid))
    plt.savefig('/data/swolk/SILVERBERG/hltau_plots/'+str(baseid)+filestem)
    #Reset for future analysis
    if not useinputbins:
        for specid in [baseidUse+1,baseidUse+2,baseidUse+3]:
            ungroup(specid)
            if binning>0:
                notice_id(specid,0.3,9.)
                mask = get_data(specid).mask
                group_counts(specid,binning,tabStops=~mask)
    return


def fit_all_obsids_hltau(baseids,tab,baselinedict,filestem,nHboundsdict,inputbins,inputbinflag=False):
    '''Fit all observations of HL Tau, taking into account whether observed with Chandra or XMM, whether observed by Chandra with or without gratings, and whether to use all three detectors on XMM.'''
    for baseid in baseids:
        startTime = time.time()
        #Set baseline assumptions--a run-of-the-mill XMM observation where MOS2 is used. Adjust from this baseline accordingly.
        a1.nH.min = 0.0
        a1.nH.max = 1000000.0
        s1.Fe.max = 1000.0
        set_analysis('energy')
        useChandra  = False
        useGratings = False
        useMOS2     = True
        
        set_model_from_dict(baselinedict)
        #Determine if fitting iron or not
        if 's1.Fe' in baselinedict.keys():
            thaw(s1.Fe)
            frozenFeUse = False
        if 's1.Fe' not in baselinedict.keys():
            freeze(s1.Fe)
            frozenFeUse = True
        print(baseid,baseid>1000000)
        
        #Evaluate whether observed by XMM or Chandra
        baseidUse = baseid*1
        testbaseid = baseid*1
        testbaseidUse = baseid*1
        if baseid > 1000000: #if XMM
            testbaseid = int(str(baseid)[-4:])
            testbaseidUse = int(str(baseid)[-4:])
        print(baseid,baseidUse,testbaseid,testbaseidUse)

        #Get ids formatted for XMM data to appropriately call all spectra.
        if testbaseid > 150:
            testbaseidUse = testbaseid*10
            baseidUse = baseidUse*10
            print(baseid,baseidUse)
            if testbaseidUse < 10000 and testbaseid<1000:
                baseidUse *= 10
                testbaseidUse *= 10
                print(baseid,baseidUse)
            if testbaseid in [1001,1101,1201,1301]:
                baseidUse = baseid*100
                print(baseid,baseidUse)
               
        if baseidUse < 1000000:
            useChandra = True
            if baseid not in [20906,18915]: #if not an observation where gratings were in place
                useGratings = True
        if (testbaseid == '401') or (str(baseid) == '109060301'): #if one of the obsids where XZ Tau bleeds into HL Tau in MOS2
            useMOS2 = False

        #Fit a no-gratings Chandra observation
        if useChandra and (not useGratings):
            print(baseid,baseidUse,(baseid in [20906,18915]))
            if 's1.Fe' in baselinedict.keys():
                s1.Fe.max = 10.
                
            #Fit and get results
            fit(baseidUse)
            fitresultsdict = get_results_dict(str(baseid),get_fit_results())
            #Check if enough degrees of freedom are available to appropriately fit iron. If not, refit with fixed iron abundance.
            if ('s1.Fe' in baselinedict.keys()) and (get_fit_results().dof < 1):
                print('Too few dof. Refitting with XMM best-fit Fe.')
                s1.Fe = baselinedict['s1.Fe']
                freeze(s1.Fe)
                frozenFeUse = True
                fit(baseidUse)#,baseidUse+1,baseidUse+2)
                fitresultsdict = get_results_dict(str(baseid),get_fit_results())
            #Get uncertainties
            conf(baseidUse)
            confresultsdict = get_confresults_dict(fitresultsdict,get_conf_results())
            
            #Re-fit with fixed iron abundance in cases where iron is free, but fitting puts iron at maximum or minimum boundaries.
            if ('s1.Fe' in baselinedict.keys()) and (not frozenFeUse) and (((confresultsdict['s1.Fe-'] is None) and (confresultsdict['s1.Fe+'] is None)) or (((confresultsdict['s1.Fe'] == 0.) or (confresultsdict['s1.Fe'] == 10.)) and (confresultsdict['s1.Fe']-np.abs(confresultsdict['s1.Fe-']) == 0.))):
                print('Poorly constrained Fe. Refitting with XMM best-fit Fe.')
                s1.Fe = baselinedict['s1.Fe']
                freeze(s1.Fe)
                frozenFeUse = True
                fit(baseidUse)
                fitresultsdict = get_results_dict(str(baseid),get_fit_results())
                conf(baseidUse)
                confresultsdict = get_confresults_dict(fitresultsdict,get_conf_results()) 

        
        if useChandra and useGratings: #if it's a gratings observation
            #Set up iron bounds if we're fitting iron
            if 's1.Fe' in baselinedict.keys():
                s1.Fe.max = 10.
                
            #Fit zeroth order, first-order MEG, first-order HEG
            fit(baseidUse,baseidUse+1,baseidUse+2)
            fitresultsdict = get_results_dict(str(baseid),get_fit_results())
            
            #Adjust and refit if 
            if ('s1.Fe' in baselinedict.keys()) and (get_fit_results().dof < 1):
                print('Too few dof. Refitting with XMM best-fit Fe.')
                s1.Fe = baselinedict['s1.Fe']
                freeze(s1.Fe)
                frozenFeUse = True
                fit(baseidUse,baseidUse+1,baseidUse+2)
                fitresultsdict = get_results_dict(str(baseid),get_fit_results())
            conf(baseidUse,baseidUse+1,baseidUse+2)
        #    confesultsdict = get_confresults_dict
            confresultsdict = get_confresults_dict(fitresultsdict,get_conf_results())
            if ('s1.Fe' in baselinedict.keys()):
                Femindif = False
                if not (confresultsdict['s1.Fe-'] is None):
                    Femindif = (confresultsdict['s1.Fe']-np.abs(confresultsdict['s1.Fe-']) == 0.)
            if ('s1.Fe' in baselinedict.keys()) and (not frozenFeUse) and (((confresultsdict['s1.Fe-'] is None) and (confresultsdict['s1.Fe+'] is None)) or (((confresultsdict['s1.Fe'] == 0.) or (confresultsdict['s1.Fe'] == 1.5)) and Femindif)):
                print('Poorly constrained Fe. Refitting with XMM best-fit Fe.')
                s1.Fe = baselinedict['s1.Fe']
                freeze(s1.Fe)
                frozenFeUse = True
                fitresultsdict = get_results_dict(str(baseid),get_fit_results())
                conf(baseidUse,baseidUse+1,baseidUse+2)
                confresultsdict = get_confresults_dict(fitresultsdict,get_conf_results()) 

        # If using all three XMM detectors
        if (not useChandra) and (useMOS2):
            fit(baseidUse+1,baseidUse+2,baseidUse+3)
            fitresultsdict = get_results_dict(str(baseid),get_fit_results())
            conf(baseidUse+1,baseidUse+2,baseidUse+3)
            confresultsdict = get_confresults_dict(fitresultsdict,get_conf_results())

        # If only using MOS1 and PN
        if (not useMOS2):           
            fit(baseidUse+1,baseidUse+3)
            fitresultsdict=get_results_dict(str(baseid),get_fit_results())
            conf(baseidUse+1,baseidUse+3)
            confresultsdict = get_confresults_dict(fitresultsdict,get_conf_results())
            
        confresultsdict = get_confresults_dict(fitresultsdict,get_conf_results())

        #Fix formatting for numbers in iron abundance if necessary
        if frozenFeUse:
            confresultsdict['s1.Fe'] = s1.Fe.val
            confresultsdict['s1.Fe-'] = 0.
            confresultsdict['s1.Fe+'] = 0.

        #Plot the data with model comparison
        if (not useChandra):
            plot_xmm_data(baseid,baseidUse,filestem,plotMOS2=useMOS2,binning=inputbins,useinputbins=inputbinflag)
        if useChandra:
            plot_chandra_data(baseid,baseidUse,filestem,gratings=useGratings,binning=inputbins,useinputbins=inputbinflag)

        #Write to table
        tab.add_row([confresultsdict[x] for x in tab.colnames])
        print(baseid,'done. Elapsed time:',time.time()-startTime)
    return tab


def get_xmm_absorbed_flux_info(obsid,noMOS2=False):
    '''Using already-set model information to evaluate source energy flux between 0.5 and 8 keV for a given XMM observation.'''
    idUse = int(obsid+'03')
    print(idUse)
    print(s1)
    absFlux = calc_energy_flux(0.5,8.,id=idUse) #Calculate the energy flux
    print(absFlux)
    
    useUncorrelated = False
    
    if noMOS2: #If only MOS1 and PN
        print(idUse-2,idUse)
        #Get a distribution of fluxes by varying the model parameters according to their uncertainties and finding flux 1000 times
        try:
            #Assume correlated uncertainties
            testDistAbs = sample_energy_flux(0.5,8.,id=idUse,otherids=[idUse-2],correlated=True,num=1000)
        #Try with uncorrelated uncertainties if correlated fails
        except TypeError:
            print('Correlated failed. Trying with uncorrelated.')
            useUncorrelated = True
            testDistAbs = sample_energy_flux(0.5,8.,id=idUse,otherids=[idUse-2],correlated=False,num=1000)
        except np.linalg.LinAlgError as e:
            print(str(e))
            print('Correlated failed. Trying with uncorrelated.')
            useUncorrelated = True
            testDistAbs = sample_energy_flux(0.5,8.,id=idUse,otherids=[idUse-2],correlated=False,num=1000)
        #Return absorbed flux, uncertainty from 16th and 84th percentile values of the sample distribution, whether to use correlated
        return absFlux,np.percentile(testDistAbs[:,0],16) - absFlux, np.percentile(testDistAbs[:,0],84) - absFlux, useUncorrelated

    #For obsids where you DO use MOS2, the same pattern
    print(idUse-2,idUse-1,idUse)
    try:
        testDistAbs = sample_energy_flux(0.5,8.,id=idUse,otherids=[idUse-2,idUse-1],correlated=True,num=1000)
    except TypeError:
        print('Correlated failed. Trying with uncorrelated.')
        useUncorrelated = True
        testDistAbs = sample_energy_flux(0.5,8.,id=idUse,otherids=[idUse-2,idUse-1],correlated=False,num=1000)
    except np.linalg.LinAlgError as e:
        print(str(e))
        print('Correlated failed. Trying with uncorrelated.')
        useUncorrelated = True
        testDistAbs = sample_energy_flux(0.5,8.,id=idUse,otherids=[idUse-2,idUse-1],correlated=False,num=1000)
    return absFlux, np.percentile(testDistAbs[:,0],16) - absFlux, np.percentile(testDistAbs[:,0],84) - absFlux, useUncorrelated

    
def get_chandra_absorbed_flux_info(obsid,useGratings=False,frozenFe=False):
    idUse = int(obsid+'0')
    if frozenFe:
        freeze(s1.Fe)

    #Get the energy flux from the best-fit parameters
    fluxAbs = calc_energy_flux(0.5,8.,id=idUse)
    useUncorrelated = False

    #If gratings, fit zeroth order, first-order MEG, first-order HEG, as in XMM
    if useGratings: 
        try:
            testDistAbs = sample_energy_flux(0.5,8.,id=idUse,otherids=[idUse+1,idUse+2],correlated=True,num=1000)
        except TypeError:
            print('Correlated failed. Trying with uncorrelated.')
            useUncorrelated = True
            testDistAbs = sample_energy_flux(0.5,8.,id=idUse,otherids=[idUse+1,idUse+2],correlated=False,num=1000)
        except np.linalg.LinAlgError as e:
            print(str(e))
            print('Correlated failed. Trying with uncorrelated.')
            useUncorrelated = True
            testDistAbs = sample_energy_flux(0.5,8.,id=idUse,otherids=[idUse+1,idUse+2],correlated=False,num=1000)            
        return(fluxAbs, np.percentile(testDistAbs[:,0],16) - fluxAbs, np.percentile(testDistAbs[:,0],84) - fluxAbs, useUncorrelated)

    #If no gratings, just use zeroth-order
    try:
        testDistAbs = sample_energy_flux(0.5,8.,id=idUse,correlated=True,num=1000)
    except TypeError: #One error that comes up along the way.
        print('Correlated failed. Trying with uncorrelated.')
        useUncorrelated = True
        testDistAbs = sample_energy_flux(0.5,8.,id=idUse,correlated=False,num=1000)
    except np.linalg.LinAlgError as e: #Another error that comes up along the way.
        print(str(e))
        print('Correlated failed. Trying with uncorrelated.')
        useUncorrelated = True
        testDistAbs = sample_energy_flux(0.5,8.,id=idUse,otherids=[idUse+1,idUse+2],correlated=False,num=1000)
    return(fluxAbs, np.percentile(testDistAbs[:,0],16) - fluxAbs, np.percentile(testDistAbs[:,0],84) - fluxAbs, useUncorrelated)


def get_xmm_unabsorbed_flux_info(obsid,modelUse,noMOS2=False,useUncorrelated=False):
    '''Get unabsorbed fluxes, based on emitting component s1.'''
    idUse = int(obsid+'03')
    print(idUse)
    print(s1)
    unabsFlux = calc_energy_flux(0.5,8.,model=(modelUse),id=idUse)
    print(unabsFlux)
    
    if noMOS2:
        print(idUse-2,idUse)
        if not useUncorrelated: #Try with correlated if that didn't fail previously
            try:
                testDistUnabs = sample_energy_flux(0.5,8.,id=idUse,model=(modelUse),otherids=[idUse-2],correlated=True,num=1000)
            except TypeError:
                print('Correlated failed. Trying with uncorrelated.')
                useUncorrelated = True
            except np.linalg.LinAlgError as e:
                print(str(e))
                print('Correlated failed. Trying with uncorrelated.')
                useUncorrelated = True

        if useUncorrelated: #Or just use uncorrelated if need be.
            print('Trying with uncorrelated.')
            testDistUnabs = sample_energy_flux(0.5,8.,id=idUse,model=(modelUse),otherids=[idUse-2],correlated=False,num=1000)
        return (unabsFlux, np.percentile(testDistUnabs[:,0],16) - unabsFlux, np.percentile(testDistUnabs[:,0],84) - unabsFlux)

    #While using MOS2
    print(idUse-2,idUse-1,idUse)
    if not useUncorrelated:
        try:
            testDistUnabs = sample_energy_flux(0.5,8.,id=idUse,model=(modelUse),otherids=[idUse-2,idUse-1],correlated=True,num=1000)
        except TypeError:
            print('Correlated failed. Trying with uncorrelated.')
            useUncorrelated = True
        except np.linalg.LinAlgError as e:
            print(str(e))
            print('Correlated failed. Trying with uncorrelated.')
            useUncorrelated = True
    if useUncorrelated:
        print('Trying with uncorrelated.')
        testDistUnabs = sample_energy_flux(0.5,8.,id=idUse,model=(modelUse),otherids=[idUse-2,idUse-1],correlated=False,num=1000)
    return (unabsFlux, np.percentile(testDistUnabs[:,0],16) - unabsFlux, np.percentile(testDistUnabs[:,0],84) - unabsFlux)

    
def get_chandra_unabsorbed_flux_info(obsid,modelUse,useGratings=False,frozenFe=False,useUncorrelated=False):
    '''Get unabsorbed flux from a Chandra observation, with flags for use of gratings spectra, and input from absorbed flux as to
    whether to default to uncorrelated variables or not.'''
    
    idUse = int(obsid+'0') #Set obsid

    if frozenFe:
        freeze(s1.Fe) #Freeze iron if there aren't enough degrees of freedom to do Fe

    print(s1)
    unabsFlux = calc_energy_flux(0.5,8.,id=idUse,model=(modelUse)) #Calculate the unabsorbed flux over the 0.5-8 keV range
    print(unabsFlux)

    #Get a sample of unabsorbed fluxes from which to estimate an uncertainty on the unabsorbed flux derived above
    #Processing for gratings spectra
    if useGratings:
        if not useUncorrelated: #Testing with correlated variables, with flags for running into failures of that method.
            try:
                testDistUnabs = sample_energy_flux(0.5,8.,id=idUse,model=(modelUse),otherids=[idUse+1,idUse+2],correlated=True,num=1000)
            except TypeError:
                print('Correlated failed. Trying with uncorrelated.')
                useUncorrelated = True
            except np.linalg.LinAlgError as e:
                print(str(e))
                print('Correlated failed. Trying with uncorrelated.')
                useUncorrelated = True
        if useUncorrelated: #Uncorrelated variables--can default to confidence intervals rather than covariance
            print('Trying with uncorrelated.')
            testDistUnabs = sample_energy_flux(0.5,8.,id=idUse,model=(modelUse),otherids=[idUse+1,idUse+2],correlated=False,num=1000)
        return(unabsFlux, np.percentile(testDistUnabs[:,0],16) - unabsFlux, np.percentile(testDistUnabs[:,0],84) - unabsFlux)

    if not useUncorrelated:
        try:
            testDistUnabs = sample_energy_flux(0.5,8.,id=idUse,model=(modelUse),correlated=True,num=1000)
        except TypeError:
            print('Correlated failed. Trying with uncorrelated.')
            useUncorrelated = True
        except np.linalg.LinAlgError as e:
            print(str(e))
            print('Correlated failed. Trying with uncorrelated.')
            useUncorrelated = True
    if useUncorrelated:
        print('Trying with uncorrelated.')
        testDistUnabs = sample_energy_flux(0.5,8.,id=idUse,model=(modelUse),correlated=False,num=1000)
    return(unabsFlux, np.percentile(testDistUnabs[:,0],16) - unabsFlux, np.percentile(testDistUnabs[:,0],84) - unabsFlux)


def get_flux_info(intab,modelUse):
    '''Run through a table and get absorbed and unabsorbed fluxes for everything in said table.'''
    checkThawedFe = False
    
    tab = copy.deepcopy(intab)
    tab['FluxAbs'] = 0.
    tab['FluxAbs'] = 0.
    tab['FluxAbsLow'] = 0.
    tab['FluxAbsHigh'] = 0.
    tab['FluxUnabs'] = 0.
    tab['FluxUnabsLow'] = 0.
    tab['FluxUnabsHigh'] = 0.

    thawedFeTable = ('s1.Fe' in tab.colnames)
    if not thawedFeTable:
        checkFrozenFe = True

    #Process the first row, which is a joint fit of the three faint obsids from 2020 that do not have XZ Tau flaring
    startTime = time.time()
    print(tab['obsid'][0])
    set_model_from_dict({col: tab[col][0] for col in tab.colnames if '.' in col and ('+' not in col) and ('-' not in col)})
    print(s1)
    tab['FluxAbs'][0] = calc_energy_flux(0.5,8.,id=86504050103)
    testDistAbs = sample_energy_flux(0.5,8.,id=86504050103,
                                     otherids=[86504030101,86504030102,86504030103,86504070101,86504070102,86504070103,
                                               86504050101,86504050102],correlated=True,num=1000)
    tab['FluxAbsLow'][0] = np.percentile(testDistAbs[:,0],16) - tab['FluxAbs'][0]
    tab['FluxAbsHigh'][0] = np.percentile(testDistAbs[:,0],84) - tab['FluxAbs'][0]
    print(tab['obsid'][0],'absorbed done. Elapsed time',time.time()-startTime)
    
    tab['FluxUnabs'][0] = calc_energy_flux(0.5,8.,model=(modelUse),id=86504050103)
    testDistUnabs = sample_energy_flux(0.5,8.,model=(modelUse),id=86504050103,
                                           otherids=[86504030101,86504030102,86504030103,86504070101,86504070102,86504070103,
                                                     86504050101,86504050102],correlated=True,num=1000)
    tab['FluxUnabsLow'][0] = np.percentile(testDistUnabs[:,0],16) - tab['FluxUnabs'][0]
    tab['FluxUnabsHigh'][0] = np.percentile(testDistUnabs[:,0],84) - tab['FluxUnabs'][0]
    print(tab['obsid'][0],'unabsorbed done. Elapsed time',time.time()-startTime)

    #Loop through the individual obsids
    for i in range(1,len(tab['obsid'])):
        startTime = time.time()
        freeze(s1.Fe)
        if thawedFeTable:
            thaw(s1.Fe)
            checkFrozenFe = False
        print(tab['obsid'][i]) #Setting the model for each obsid
        set_model_from_dict({col: tab[col][i] for col in tab.colnames if '.' in col and ('+' not in col) and ('-' not in col)})
        pprint.pprint({col: tab[col][i] for col in tab.colnames if '.' in col and ('+' not in col) and ('-' not in col)})

        if thawedFeTable and ((tab['s1.Fe-'][i] == 0.) and (tab['s1.Fe+'][i] == 0.)):
            checkFrozenFe = True

        #Processing for if it's an observation with the MOS2 PSF from XZ Tau bleeding directly into the HL Tau PSF.
        if int(tab['obsid'][i]) in [865040401,109060301]:
            print(a1)
            print(s1)
            tab['FluxAbs'][i], tab['FluxAbsLow'][i], tab['FluxAbsHigh'][i], uncorXMM = get_xmm_absorbed_flux_info(tab['obsid'][i],noMOS2=True)
            print(tab['obsid'][i],'absorbed done. Elapsed time',time.time()-startTime)

            print(a1)
            print(s1)
            tab['FluxUnabs'][i], tab['FluxUnabsLow'][i], tab['FluxUnabsHigh'][i] = get_xmm_unabsorbed_flux_info(tab['obsid'][i],modelUse,noMOS2=True, useUncorrelated=uncorXMM)
            print(tab['obsid'][i],'done. Elapsed time',time.time()-startTime)
            continue

        #Processing for other XMM observations
        if int(tab['obsid'][i]) > 1000000:
            tab['FluxAbs'][i], tab['FluxAbsLow'][i], tab['FluxAbsHigh'][i], uncorXMM = get_xmm_absorbed_flux_info(tab['obsid'][i])
            print(tab['obsid'][i],'absorbed done. Elapsed time',time.time()-startTime)
            
            tab['FluxUnabs'][i], tab['FluxUnabsLow'][i], tab['FluxUnabsHigh'][i] = get_xmm_unabsorbed_flux_info(tab['obsid'][i],modelUse,useUncorrelated=uncorXMM)
            print(tab['obsid'][i],'done. Elapsed time',time.time()-startTime)
            continue

        #Processing Chandra observations without gratings
        if (int(tab['obsid'][i]) in [20906,18915]):
            tab['FluxAbs'][i], tab['FluxAbsLow'][i], tab['FluxAbsHigh'][i], uncorChandra = get_chandra_absorbed_flux_info(tab['obsid'][i],frozenFe=checkFrozenFe)
            print(tab['obsid'][i],'absorbed done. Elapsed time',time.time()-startTime)
            
            tab['FluxUnabs'][i], tab['FluxUnabsLow'][i], tab['FluxUnabsHigh'][i] = get_chandra_unabsorbed_flux_info(tab['obsid'][i],modelUse,frozenFe=checkFrozenFe,useUncorrelated=uncorChandra)
            print(tab['obsid'][i],'done. Elapsed time',time.time()-startTime)
            continue

        #Processing Chandra observations with gratings
        if (int(tab['obsid'][i]) > 10000) and (int(tab['obsid'][i]) < 1000000):
            tab['FluxAbs'][i], tab['FluxAbsLow'][i], tab['FluxAbsHigh'][i], uncorChandra = get_chandra_absorbed_flux_info(tab['obsid'][i], useGratings=True, frozenFe=checkFrozenFe)
            print(tab['obsid'][i],'absorbed done. Elapsed time',time.time()-startTime)
            
            tab['FluxUnabs'][i], tab['FluxUnabsLow'][i], tab['FluxUnabsHigh'][i] = get_chandra_unabsorbed_flux_info(tab['obsid'][i], modelUse, useGratings=True, frozenFe=checkFrozenFe,useUncorrelated=uncorChandra)
            print(tab['obsid'][i],'done. Elapsed time',time.time()-startTime) 

    return tab