AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = libjetanal.la

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(EICDIRECTORY)/include \
  -I$(ROOTSYS)/include \
  -I`root-config --incdir`

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(EICDIRECTORY)/lib \
  -L$(ROOTSYS)/lib \
  `root-config --libs` \
  -L/afs/rhic.bnl.gov/eic/restructured/env/pro/lib

ACLOCAL_AMFLAGS=-I m4

CXX=g++
CXXFLAGS = -g -Wall -std=c++11

pkginclude_HEADERS = \
  JetAnalysis.h \
  EventShape.h

if MAKEROOT6
  pcmdir = $(libdir)
  nobase_dist_pcm_DATA = \
    EventShape_Dict_rdict.pcm \
    JetAnalysis_Dict_rdict.pcm 
else
  ROOT5_IO_DICTS = \
    JetAnalysis_Dict.cc \
    EventShape_Dict.cc
endif

libjetanal_la_SOURCES = \
  $(ROOT5_IO_DICTS) \
  EventShape.cc \
  JetAnalysis.cc 

libjetanal_la_LIBADD = \
  -lfastjet \
  -lfastjettools \
  -leicsmear

# Rule for generating CINT dictionaries from class headers.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@  -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests
#
# BUILT_SOURCES = testexternals.cc
#
# noinst_PROGRAMS = \
#   testexternals_track_io
#
#   #testexternals_track_util_SOURCES = testexternals.cc
#testexternals_track_util_LDADD = libtrack_util.la

testexternals_jetanal_SOURCES = testexternals.cc
testexternals_jetanal_LDADD = libjetanal.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

################################################



clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
