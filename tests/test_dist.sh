#! /bin/sh
#_______________________________________________________________________________
# test OPTICS on dist data

#../src/optics_dist --datafile dist.dat || exit 1
valgrind --leak-check=full --show-leak-kinds=all ../src/optics_dist --datafile dist.dat || exit 1
