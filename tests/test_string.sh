#! /bin/sh
#_______________________________________________________________________________
# test OPTICS on string data

../src/optics_str --datafile string.dat --minpts 50 --eps 50 || exit 1
