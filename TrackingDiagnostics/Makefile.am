AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

pkginclude_HEADERS = \
  DistortionAnalysis.h \
  PHTrackPruner.h

lib_LTLIBRARIES = \
  libDistortionAnalysis.la

libDistortionAnalysis_la_SOURCES = \
  DistortionAnalysis.cc \
  PHTrackPruner.cc

libDistortionAnalysis_la_LIBADD = \
  -lphool \
  -lfun4all \
  -lg4eval \
  -lg4eval_io \
  -lphg4hit \
  -lg4detectors_io \
  -lmvtx \
  -lmvtx_io \
  -lintt \
  -lintt_io \
  -lmicromegas_io \
  -lSubsysReco \
  -ltpc_io \
  -ltrack_io \
  -ltrackeralign \
  -ltrackbase_historic_io

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libDistortionAnalysis.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
