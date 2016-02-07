#! /bin/sh
#_______________________________________________________________________________
#  test OPTICS on xyz coordinate file

../src/optics_xyz --datafile point.dat --minpts 2 || exit 1
