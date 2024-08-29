from ciao_contrib.runtool import *
from create_filtered_evtfile import create_filtered_evtfile

def extract_hltau_lightcurve(obsid, chip,binsize):
    #Create eventfile filtered on 0.5-8 keV
    create_filtered_evtfile(obsid, '/data/swolk/SILVERBERG/hltau_data/')
    
    dmextract.punlearn()
    dmextract.infile  = '/data/swolk/SILVERBERG/hltau_data/'+str(obsid)+'/repro/acisf'+str(obsid)+'_repro_evt2_500-8000.fits[ccd_id='+str(chip)+',sky=region(/data/swolk/SILVERBERG/hltau_xztau_scripts/hltausrc.reg)][bin time=::'+str(binsize)+']'
    dmextract.outfile = '/data/swolk/SILVERBERG/hltau_data/'+str(obsid)+'/repro/hltau_sub_lc_500_8000.fits'
    dmextract.bkg     = '/data/swolk/SILVERBERG/hltau_data/'+str(obsid)+'/repro/acisf'+str(obsid)+'_repro_evt2_500-8000.fits[ccd_id='+str(chip)+',sky=region(/data/swolk/SILVERBERG/hltau_xztau_scripts/xztaubkg.reg)]'
    dmextract.opt     = 'ltc1'
    dmextract.verbose = 2
    dmextract.clobber = True
    
    a = dmextract()
    print(a)
    return
