#SUBDIRS = 0865040401 0865040601 
#SUBDIRS = 0865040401 0865040601 0865040501 0865040701 0865040201 0865040301   `
#SUBDIRS = 0200810901 0094810701 0200811001 0200810201 0094810201 0200810501 0200811301
#SUBDIRS = 0865040301

#SUBDIRS = 0865040401 0865040601 0865040501 0865040701 0865040201 0865040301   `
#SUBDIRS = 0865040401 0865040601 0865040501 0865040701 0865040201 0865040301 0200810401 0200811101 0200810501 0200811201 0200810601 0109060301 0200810701 0109060501 0200810801 0200810201 0200810901 0200810301 0200811001 
#SUBDIRS = 0762360101
#SUBDIRS = 0200811201

SUBDIRS = 0109060301 #0200810201 0200810301 0200810401 0200810501 0200810601 0200810701 0200810801 0200810901 0200811001 0200811101 0200811201 0200811301

ORIGIN = /data/swolk/SILVERBERG/hltau_xztau_scripts/
SRCFILES = /data/swolk/SILVERBERG/hltau_xztau_scripts/hltau_src/
FLARE = 

.PHONY: subdirs $(SUBDIRS) noflare

subdirs: $(SUBDIRS) noflare

makefile : $(ORIGIN)/makefile
	cp $< $@

defaults.mk : $(ORIGIN)/defaults.mk
	cp $< $@

$(SUBDIRS):
	# Auto update makefile in subdir from master make
	$(MAKE) -C $@ makefile -f $(ORIGIN)/master.mk
	$(MAKE) -C $@ defaults.mk -f $(ORIGIN)/master.mk
	$(MAKE) -C $@ SRCFILE=$@.src EPIC_images
	cp $(SRCFILES)/$@.src $@/
	#$(MAKE) -C $@ SRCFILE=$@.src EPIC_lc
	$(MAKE) -C $@ clean_EPIC
	$(MAKE) -C $@ SRCFILE=$@.src EPIC

# deal with those ObsIDs which need anoflare filtering and have special noflare src files
# in this setup it only works if that is only one dir!
#otherwise need some string procession to get the _noflare part removed from dir name
noflare : $(FLARE)
	#cp $(SRCFILES)/$<_noflare.src $</
	#$(MAKE) -C $< SRCFILE=$<_noflare.src EPIC

# in /data/swolk/SILVERBERG/hltau_data/
# make subdirs -f /data/swolk/SILVERBERG/hltau_xztau_scripts/master.mk

# in /media/MAX/moritz/obs/XMM/
# make subdirs -f /nfs/melkor/d1/guenther/soft/XMM-data-reduction-scripts/master.mk -K -j 3
# -j 3: Do parallel, but max 3 processes at a time. I have dual core, if I use much more, than it's swapping all the time.bstracts
#  0200810201  0200810901  0865040401  package_24241_210628145418.tar.0
#  0200810301  0200811001  0865040501
