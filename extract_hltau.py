from ciao_contrib.runtool import *
from get_skycoords import get_skycoords
import os
import subprocess
import glob

def extract_hltau(srcid,ra,dec):
    '''Extract TG spectrum for HL Tau from a given obs id. Inputs: srcid, ra, dec'''
    srcstring = str(srcid)
    
    print('Now running',srcstring)
    print('Running dmcoords')    
    #Set up dmcoords to get sky coordinates of source from RA and DEC
    dmcoords.punlearn()
    dmcoords.option = 'cel'
    dmcoords.ra = ra
    dmcoords.dec = dec
    
    #Get skycoords from RA and DEC
    testdmcoords = dmcoords(srcstring+'/repro/acisf'+srcstring+'_repro_evt2.fits',verbose=2)
    skyx, skyy = get_skycoords(testdmcoords)
    
    #Hard-code file names for future use
    #lsString = "ls "+srcstring+"/secondary/*.fits"
    #testfilenames = subprocess.run(["ls",srcstring+"/secondary/*.fits"],shell=True)
    #print(testfilenames)

    #proc = subprocess.Popen('ls '+srcstring+'/secondary/*.fits',shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    #testfilenames,err = proc.communicate()
    #print(testfilenames)

    testfilenames = glob.glob(srcstring+'/secondary/*.fits')
    print(testfilenames)

    evt1filename        = [x for x in testfilenames if x[-9:] == 'evt1.fits'][0]
    tgdetectoutfilename = srcstring+'/acis_'+srcstring+'_hltau.fits'
    maskfilename        = srcstring+'/acis_'+srcstring+'_hltau_evt1_L1a.fits'
    evt1afilename       = srcstring+'/acis_'+srcstring+'_hltau_evt1a.fits'
    filterfile1name     = srcstring+'/acis_'+srcstring+'_hltau_flt1_evt1a.fits'
    evt2filename        = srcstring+'/acis_'+srcstring+'_hltau_evt2.fits'
    filter2filename     = [x for x in testfilenames if x[-9:] == 'flt1.fits'][0]
    
    print('Running tgdetect2')
    #Use those coordinates to set up tgdetect2, then run tgdetect2
    tgdetect2.punlearn()
    tgdetect2.zo_pos_x = skyx
    tgdetect2.zo_pos_y = skyy
    tgdetect2.clobber = True
    #if tgdetect2.clobber != clobber:
        
    tgdetect2.infile = evt1filename
    tgdetect2.outfile = tgdetectoutfilename
    tgdetect2.verbose = 2
    
    a = tgdetect2()
    print(a)
    
    print('Creating mask')
    tg_create_mask.punlearn()
    tg_create_mask.infile  = evt1filename
    tg_create_mask.outfile = maskfilename
    tg_create_mask.input_pos_tab = tgdetectoutfilename
    tg_create_mask.verbose = 2
    tg_create_mask.clobber = True
    
    b = tg_create_mask()
    print(b)
    
    print('Resolving events')
    #tg_resolve_events
    tg_resolve_events.punlearn()
    tg_resolve_events.infile = evt1filename
    tg_resolve_events.outfile = evt1afilename
    tg_resolve_events.regionfile = maskfilename
    tg_resolve_events.acaofffile = srcstring+'/repro/pcadf'+srcstring+'_000N001_asol1.fits'
    tg_resolve_events.verbose = 2
    tg_resolve_events.clobber = True
    
    c = tg_resolve_events()
    print(c)
    
    print('Filtering events')
    #Filter events, first for grade and status
    dmcopy.punlearn()
    dmcopy.infile = evt1afilename+'[EVENTS][grade=0,2,3,4,6,status=0]'
    dmcopy.outfile = filterfile1name
    dmcopy.verbose = 2
    dmcopy.clobber = True
    
    d = dmcopy()
    print(d)
    
    dmappend.punlearn()
    dmappend.infile = evt1afilename+'[region][subspace -time]'
    dmappend.outfile = filterfile1name
    
    d1 = dmappend()
    print(d1)
    
    print('Second filter')
    #Second filter
    dmcopy.punlearn()
    #dmcopy.infile = filterfile1name+'[EVENTS][@'+srcstring+'/secondary/acisf'+srcstring+'_000N002_flt1.fits][cols -phas]'
    dmcopy.infile  = filterfile1name+'[EVENTS][@'+filter2filename+'][cols -phas]'
    dmcopy.outfile = evt2filename
    dmcopy.verbose = 2
    dmcopy.clobber = True    

    d2 = dmcopy()
    print(d2)
    
    dmappend.punlearn()
    
    dmappend.infile=evt1afilename+'[region][subspace -time]'
    dmappend.outfile = evt2filename
    dmappend.verbose = 2
    
    d3 = dmappend()
    print(d3)
    
    print('Running tgextract')
    #tgextract
    tgextract.punlearn()
    tgextract.infile = evt2filename
    tgextract.outfile = srcstring+'/repro/acis_'+srcstring+'_hltau_pha2.fits'
    tgextract.verbose = 2
    tgextract.clobber = True   
 
    f = tgextract()
    print(f)

    
    #Make ARF and RMF for observation
    mktgresp.punlearn()
    mktgresp.infile = srcstring+'/repro/acis_'+srcstring+'_hltau_pha2.fits'
    mktgresp.evtfile = evt2filename
    mktgresp.outroot = srcstring+'/repro/tg_hltau/acis_'+srcstring+'_hltau'
    mktgresp.verbose = 2
    mktgresp.clobber = True

    g = mktgresp()
    print(g)
    print('\nDone')
    return
