##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include \
  -isystem$(OPT_SPHENIX)/include


AM_LDFLAGS = \
  -L$(libdir) \
  -L$(ROOTSYS)/lib \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

# List of shared libraries to produce
lib_LTLIBRARIES = \
  libtpccalib_io.la \
  libtpccalib.la

libtpccalib_io_la_LIBADD = \
  -lphool

libtpccalib_la_LIBADD = \
  libtpccalib_io.la \
  -lFROG \
  -lffaobjects \
  -lphool \
  -lodbc++ \
  -lSubsysReco \
  -lg4detectors_io \
  -ltrack_io \
  -ltrackbase_historic_io \
  -ltrack_reco \
  -ltpc_io

pkginclude_HEADERS = \
  TpcDirectLaserReconstruction.h \
  TpcSpaceChargeMatrixContainer.h \
  TpcSpaceChargeMatrixContainerv1.h \
  TpcSpaceChargeMatrixContainerv2.h \
  TpcSpaceChargeMatrixContainer1D.h \
  TpcSpaceChargeMatrixContainer2D.h \
  TpcSpaceChargeMatrixInversion.h \
  TpcSpaceChargeMatrixInversion1D.h \
  TpcSpaceChargeMatrixInversion2D.h \
  TpcSpaceChargeReconstruction.h \
  TpcSpaceChargeReconstructionHelper.h \
  PHTpcResiduals.h \
  TpcCentralMembraneMatching.h \
  TpcLaminationFitting.h


ROOTDICTS = \
  TpcSpaceChargeMatrixContainer_Dict.cc \
  TpcSpaceChargeMatrixContainerv1_Dict.cc \
  TpcSpaceChargeMatrixContainerv2_Dict.cc \
  TpcSpaceChargeMatrixContainer1D_Dict.cc \
  TpcSpaceChargeMatrixContainer2D_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  TpcSpaceChargeMatrixContainer_Dict_rdict.pcm \
  TpcSpaceChargeMatrixContainerv1_Dict_rdict.pcm \
  TpcSpaceChargeMatrixContainerv2_Dict_rdict.pcm \
  TpcSpaceChargeMatrixContainer1D_Dict_rdict.pcm \
  TpcSpaceChargeMatrixContainer2D_Dict_rdict.pcm

libtpccalib_io_la_SOURCES = \
  $(ROOTDICTS) \
  TpcSpaceChargeMatrixContainerv1.cc \
  TpcSpaceChargeMatrixContainerv2.cc \
  TpcSpaceChargeMatrixContainer1D.cc \
  TpcSpaceChargeMatrixContainer2D.cc

libtpccalib_la_SOURCES = \
  TpcDirectLaserReconstruction.cc \
  TpcSpaceChargeMatrixInversion.cc \
  TpcSpaceChargeMatrixInversion1D.cc \
  TpcSpaceChargeMatrixInversion2D.cc \
  TpcSpaceChargeReconstruction.cc \
  TpcSpaceChargeReconstructionHelper.cc \
  PHTpcResiduals.cc \
  TpcCentralMembraneMatching.cc \
  TpcLaminationFitting.cc

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

BUILT_SOURCES = \
  testexternals.cc

noinst_PROGRAMS = \
  testexternals_tpccalib

testexternals_tpccalib_SOURCES = testexternals.cc
testexternals_tpccalib_LDADD = libtpccalib.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

##############################################
# please add new classes in alphabetical order

clean-local:
	rm -f $(BUILT_SOURCES)
