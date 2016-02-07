#! /bin/csh
#_______________________________________________________________________________
# test OPTICS MPI on xyz coordinate file

mpirun -np 1 ../src/optics_xyz_mpi --datafile point.dat || exit 1
