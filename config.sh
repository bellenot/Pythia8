#!/bin/sh
#
# Find platform.
#
export ARCH="`uname`"
export theGcc1=`g++ --version | awk '{print$3}'`
export theGcc=`echo ${theGcc1} | awk -F . '{print $1}'`
if [ ${theGcc} = 4 ]; then
  export ARCH=${ARCH}-gcc4
fi
echo Platform is $ARCH
#
# Environment variables managing the usage of Pythia6 library. If the
# variables below are not set, the local Pythia6 library will be built
# and used.
# The default values here correspond to CERN AFS (but you may want to change
# the version and/or platform).
#
#export PYTHIA6LOCATION=/afs/cern.ch/sw/lcg/app/releases/GENSER/GENSER_1_3_0/slc3_ia32_gcc323/lib/archive
#export PYTHIA6LIBNAME="-lpythia6_326 -ldummypythia6_326 -lpdfdummypythia6_326"
#
# Environment variables for building HepMC interface library. Comment out ALL
# of them if this interface is not needed or if you don't have both CLHEP
# and HepMC in your computer. Note that this library is used by the example
# main11.
# The default values here correspond to CERN AFS (but you may want to change
# the version and/or platform).
#
#export HEPMCLOCATION=/afs/cern.ch/sw/lcg/external/HepMC/1.26/slc3_ia32_gcc323
#export CLHEPLOCATION=/afs/cern.ch/sw/lcg/external/clhep/1.9.2.2/slc3_ia32_gcc323
#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HEPMCLOCATION}/lib:${CLHEPLOCATION}/lib
