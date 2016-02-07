#! /bin/sh
#_______________________________________________________________________________
# test OPTICS on nested xyz coordinates

../src/optics_xyz --datafile nested.dat --minpts 2 || exit 1
