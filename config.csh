#!/bin/csh
#
# Find platform.
#
setenv ARCH "`uname`"
set theGcc1=`g++ --version | awk '{print$3}'`
set theGcc=`echo $theGcc1 | awk -F . '{print $1}'`
if( $theGcc == 4 ) then
  if( -e /usr/bin/gfortran ) then
    setenv ARCH ${ARCH}-gcc4
  else
    echo ERROR: compiler is gcc4 but gfortran is not installed
  endif
endif
echo Platform is $ARCH
#
# Environment variables managing the usage of Pythia6 library. If the
# variables below are not set, the local Pythia6 library will be built
# and used.
# The default values here correspond to CERN AFS (but you may want to change
# the version and/or platform).
#
#setenv PYTHIA6LOCATION /afs/cern.ch/sw/lcg/app/releases/GENSER/GENSER_1_3_0/slc3_ia32_gcc323/lib/archive
#setenv PYTHIA6LIBNAME "-lpythia6_326 -ldummy_pythia6_326 -lpdfdummy_pythia6_326"
#
# Environment variables for building HepMC interface library. Comment out ALL
# of them if this interface is not needed or if you don't have both CLHEP
# and HepMC in your computer. Note that this library is used by the example
# main11.
# The default values here correspond to CERN AFS (but you may want to change
# the version and/or platform).
#
#setenv HEPMCLOCATION /afs/cern.ch/sw/lcg/external/HepMC/1.26/slc3_ia32_gcc323
#setenv CLHEPLOCATION /afs/cern.ch/sw/lcg/external/clhep/1.9.2.2/slc3_ia32_gcc323
#if( $?LD_LIBRARY_PATH ) then
#  setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${HEPMCLOCATION}/lib:${CLHEPLOCATION}/lib
#else
#  setenv LD_LIBRARY_PATH ${HEPMCLOCATION}/lib:${CLHEPLOCATION}/lib
#endif
