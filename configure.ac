#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT(OPTICS, 1.2.0, jens@jkleinj.eu)
AC_CONFIG_SRCDIR([src/arg.c])
AC_CONFIG_HEADERS([src/config.h])

AM_INIT_AUTOMAKE

AC_SUBST(INTI_CFLAGS)
AC_SUBST(INTI_LIBS)

# Checks for programs.
AC_PROG_CC
AC_PROG_LN_S

# Checks for libraries.
AC_CHECK_LIB([m], [cos])

# Checks for header files.
AC_CHECK_HEADERS([float.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC

# Enable the use of angle coordinates
# allow Makefile.am to use DATAANG in conditional compilation and disable option by default
AM_CONDITIONAL(DATAANG, false)
# 'configure --enable-data-ang' enables DATAANG and defines it in config.h
AC_ARG_ENABLE(data-ang,
 [  --enable-data-ang  Enable data-ang],
 [AM_CONDITIONAL(DATAANG, true) AC_DEFINE(DATAANG, [], [input data in angle format])]
)

# Enable the use of strings
# allow Makefile.am to use DATASTR in conditional compilation and disable option by default
AM_CONDITIONAL(DATASTR, false)
# 'configure --enable-data-str' enables DATASTR and defines it in config.h
AC_ARG_ENABLE(data-str,
 [  --enable-data-str  Enable data-str],
 [AM_CONDITIONAL(DATASTR, true) AC_DEFINE(DATASTR, [], [input data in str format])]
)

# Enable the use of vec coordinates
# allow Makefile.am to use DATAVEC in conditional compilation and disable by default
AM_CONDITIONAL(DATAVEC, false)
# 'configure --enable-data-vec' enables DATAVEC and defines it in config.h
AC_ARG_ENABLE(data-vec,
 [  --enable-data-vec  Enable data-vec],
 [AM_CONDITIONAL(DATAVEC, true) AC_DEFINE(DATAVEC, [], [input data in vec format])]
)

# Enable the use of xyz coordinates
# allow Makefile.am to use DATAXYZ in conditional compilation and disable by default
AM_CONDITIONAL(DATAXYZ, false)
# 'configure --enable-data-xyz' enables DATAXYZ and defines it in config.h
AC_ARG_ENABLE(data-xyz,
 [  --enable-data-xyz  Enable data-xyz],
 [AM_CONDITIONAL(DATAXYZ, true) AC_DEFINE(DATAXYZ, [], [input data in xyz format])]
)

# Enable the use of distance matrix
# allow Makefile.am to use DATADIST in conditional compilation and disable by default
AM_CONDITIONAL(DATADIST, false)
# 'configure --enable-data-dist' enables DATADIST and defines it in config.h
AC_ARG_ENABLE(data-dist,
 [  --enable-data-dist  Enable data-dist],
 [AM_CONDITIONAL(DATADIST, true) AC_DEFINE(DATADIST, [], [input data in distance format])]
)

# Enable pattern matching library 
AC_ARG_WITH(pcre,
	[  --with-pcre=<path>      prefix of pcre installation (eg /usr/local)],
	[
		CPPFLAGS="$CPPFLAGS -I $withval/include"
		LDFLAGS="$LDFLAGS -L $withval/lib"
	]
)
# look for pcre pattern matching library
AC_CHECK_HEADERS(pcre.h,
	AC_SEARCH_LIBS([pcre_compile], [pcre]),
	AC_MSG_WARN([*** pcre.h not found -- consider using --with-pcre])
)

# enable debug
AC_ARG_ENABLE(debug,
 [  --enable-debug  Turn on debugging], [ac_use_debug="yes"], [ac_use_debug="no"])
if test "$ac_use_debug" = "yes"; then
    CFLAGS="-ggdb"
    AC_DEFINE(DEBUG,1,[Define to compile with DEBUG settings.])
fi

# enable profiling 
AC_ARG_ENABLE(profiling,
 [  --enable-profiling  Turn on profiling], [ac_use_profiling="yes"], [ac_use_profiling="no"])
if test "$ac_use_profiling" = "yes"; then
    CFLAGS="-pg -O2"
    AC_DEFINE(PROFILING,1,[Define to compile with PROFILING settings.])
fi

AC_CHECK_FUNCS([pow sqrt])

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 tests/Makefile])

AC_OUTPUT

