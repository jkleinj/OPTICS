#===============================================================================
# Makefile.am : for automake
# (C) 2008-2012 Jens Kleinjung and Alessandro Pandini
# Read the COPYING file for license information.
#===============================================================================

AM_CPPFLAGS = $(INTI_CFLAGS)
AM_CFLAGS = -Wall

#______________________________________________________________________________
# build for angle coordinates
if DATAANG
bin_PROGRAMS = optics_ang

optics_ang_SOURCES	= \
arg.c arg.h config.h coords_ang.c coords_ang.h \
optics.c optics.h \
safe.c safe.h sort.c sort.h vector.c vector.h

optics_ang_LDADD = $(INTI_LIBS)
endif

#______________________________________________________________________________
# build for string coordinates
if DATASTR
bin_PROGRAMS = optics_str

optics_str_SOURCES	= \
arg.c arg.h config.h coords_str.c coords_str.h \
optics.c optics.h \
safe.c safe.h sort.c sort.h vector.c vector.h

optics_str_LDADD = $(INTI_LIBS)
endif

#______________________________________________________________________________
# build for vector (arbitrary length) coordinates
if DATAVEC
bin_PROGRAMS = optics_vec

optics_vec_SOURCES	= \
arg.c arg.h config.h coords_vec.c coords_vec.h \
optics.c optics.h \
safe.c safe.h sort.c sort.h vector.c vector.h

optics_vec_LDADD = $(INTI_LIBS)
endif

#______________________________________________________________________________
# build for xyz coordinates
if DATAXYZ
bin_PROGRAMS = optics_xyz

optics_xyz_SOURCES	= \
arg.c arg.h config.h coords_xyz.c coords_xyz.h \
optics.c optics.h \
safe.c safe.h sort.c sort.h vector.c vector.h

optics_xyz_LDADD = $(INTI_LIBS)
endif

#______________________________________________________________________________
# build for distance matrix
if DATADIST
bin_PROGRAMS = optics_dist

optics_dist_SOURCES	= \
arg.c arg.h config.h coords_dist.c coords_dist.h matrix.c matrix.h \
optics.c optics.h \
safe.c safe.h sort.c sort.h vector.c vector.h

optics_dist_LDADD = $(INTI_LIBS)
endif

#______________________________________________________________________________
EXTRA_DIST = doxygen.cfg

CLEANFILES = $(TARGET) *.o 
DISTCLEANFILES = libtool config.cache config.log

