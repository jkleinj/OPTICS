#! /bin/sh
#_______________________________________________________________________________
# test OPTICS on angle coordinate file

../src/optics_ang --datafile ang.dat --eps 100 --minpts 100 || exit 1
