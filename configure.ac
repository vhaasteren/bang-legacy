dnl ##  configure.ac -- GNU autoconf source for BAPTise, a Bayesian Analysis
dnl ##  program for Pulsar Timing Array projects
dnl ##  Copyright (c) 2007-2008 Rutger van Haasteren, Leiden Observatory
dnl ##
dnl ##  Process this file with autoconf to produce a configure script.

dnl ##  This program is free software; you can redistribute it and/or modify
dnl ##  it under the terms of the GNU General Public License as published by
dnl ##  the Free Software Foundation; either version 3, or (at your option)
dnl ##  any later version.

dnl ##  This program is distributed in the hope that it will be useful,
dnl ##  but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl ##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl ##  GNU General Public License for more details.

dnl ##  You should have received a copy of the GNU General Public License
dnl ##  along with this program; if not, write to the Free Software Foundation,
dnl ##  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.  */


dnl --------------------------------
dnl -----------------------------------------------
dnl Check for compile options
dnl -----------------------------------------------


dnl Initialization macros.
dnl --------------------------------

AC_INIT(src/banginterface.cpp)
AC_CONFIG_SRCDIR(src/banginterface.cpp)
AC_CONFIG_HEADERS(config.h)
AC_PREREQ(2.57)


dnl -----------------------------------------------
dnl Package name and version number (user defined)
dnl -----------------------------------------------

GENERIC_LIBRARY_NAME=libbang
GENERIC_PACKAGE_NAME=bang

#release versioning
GENERIC_MAJOR_VERSION=0
GENERIC_MINOR_VERSION=1
GENERIC_MICRO_VERSION=0

#API version (often = GENERIC_MAJOR_VERSION.GENERIC_MINOR_VERSION)
GENERIC_API_VERSION=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION
AC_SUBST(GENERIC_API_VERSION)

#shared library versioning
GENERIC_LIBRARY_VERSION=$GENERIC_MAJOR_VERSION:$GENERIC_MINOR_VERSION:$GENERIC_MICRO_VERSION
#  The usual scheme GLV=1:3:0
#                       | | |
#                +------+ | +---+
#                |        |     |
#             current:revision:age
#                |        |     |
#                |        |     +- increment if interfaces have been added
#                |        |        set to zero if interfaces have been removed
#                |        |        or changed
#                |        +- increment if source code has changed
#                |           set to zero if current is incremented
#                +- increment if interfaces have been added, removed or changed


dnl --------------------------------
dnl Package name and version number
dnl --------------------------------

AC_SUBST(GENERIC_LIBRARY_VERSION)

PACKAGE=$GENERIC_PACKAGE_NAME
AC_SUBST(GENERIC_LIBRARY_NAME)
AC_SUBST(GENERIC_PACKAGE_NAME)

GENERIC_VERSION_L=$GENERIC_MAJOR_VERSION_$GENERIC_MINOR_VERSION_$GENERIC_MICRO_VERSION
GENERIC_VERSION=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION.$GENERIC_MICRO_VERSION
GENERIC_RELEASE=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION
AC_SUBST(GENERIC_RELEASE)
AC_SUBST(GENERIC_VERSION)

VERSION=$GENERIC_RELEASE
AC_SUBST(VERSION)

AM_INIT_AUTOMAKE($PACKAGE, $VERSION)




dnl -----------------------------------------------
dnl Checks for programs.
dnl -----------------------------------------------

AC_PROG_CXX([g++])
#AC_PROG_CXX([g++34])
#AC_PROG_CXX([mpicxx])
AC_STDC_HEADERS
AC_PROG_INSTALL

AM_SANITY_CHECK
AC_LANG_CPLUSPLUS


dnl -----------------------------------------------
dnl Declare that we are building a library
dnl -----------------------------------------------

dnl AC_DISABLE_STATIC
dnl AC_DISABLE_SHARED
AC_PROG_LIBTOOL


dnl -----------------------------------------------
dnl Path to libraries that are not installed
dnl -----------------------------------------------

#CUSTOM_LIB_PATH="lib"
CUSTOM_INCLUDE_PATH="$TEMPO2/include"

CPPFLAGS="$CPPFLAGS -I$CUSTOM_INCLUDE_PATH"
#LDFLAGS="$LDFLAGS -L../$CUSTOM_LIB_PATH  -L../../$CUSTOM_LIB_PATH"

#LIBRARY_PATH="$LIBRARY_PATH:./$CUSTOM_LIB_PATH"
#LD_LIBRARY_PATH="$LD_LIBRARY_PATH:./$CUSTOM_LIB_PATH"
C_INCLUDE_PATH="$C_INCLUDE_PATH:$CUSTOM_INCLUDE_PATH"
#export LIBRARY_PATH LD_LIBRARY_PATH C_INCLUDE_PATH
export C_INCLUDE_PATH


dnl -----------------------------------------------
dnl Check for compile options
dnl -----------------------------------------------

# We need a BLAS/cblas library. These 3 are supported:
# 'gsl' is the gnu scientific library, available pre-packaged for most platforms
# 'atlas' is the automatically tunes blas library, also usually pre-packaged
# 'mkl' is the Intel MKL blas/lapack library.
AC_ARG_WITH(blaslib,
	[  --with-blaslib[=ARG]      Specify the blas library (gsl,veclib,atlas,mkl,custom default:ARG=gsl)],
	[
		if   test "x$with_blaslib" = "xgsl" ; then
			blas_lib=gsl
		elif test "x$with_blaslib" = "xatlas" ; then
			blas_lib=atlas
		elif test "x$with_blaslib" = "cblas" ; then
			blas_lib=cblas
		elif test "x$with_blaslib" = "xmkl" ; then
			blas_lib=mkl
		elif test "x$with_blaslib" = "xveclib" ; then
		  	blas_lib=veclib
		else
			blas_lib="$with_blaslib"
		fi
	],
	[
		blas_lib=gsl
	])

# We also need a lapack/linal library. These 3 are supported:
# 'gsl' is the gnu scientific library, available pre-packaged for most platforms
# 'mkl' is the Intel MKL blas/lapack library.
# 'netlib' (www.netlib.org) provides a f2c'd version of the FORTRAN lapack
#          library. Use it with the blas wrapper libraries
#
# Note that also 'clapack' exists, but this library is not complete (no
# eigenfunctions for instance). Therefore only the complete FORTRAN libraries
# above are supported.
AC_ARG_WITH(lapacklib,
	[  --with-lapacklib[=ARG]    Specify the lapack library (gsl,veclib,netlib,mkl default:ARG=lapack)],
	[
		if   test "x$with_lapacklib" = "xgsl" ; then
			lapack_lib=gsl
		elif test "x$with_lapacklib" = "xclapack" ; then
			lapack_lib=clapack
		elif test "x$with_lapacklib" = "xnetlib" ; then
			lapack_lib=netlib
		elif test "x$with_lapacklib" = "xmkl" ; then
			lapack_lib=mkl
		elif test "x$with_lapacklib" = "xlapack" ; then
			lapack_lib=lapack
		elif test "x$with_lapacklib" = "xveclib" ; then
		  	lapack_lib=veclib
		else
			lapack_lib="$with_lapacklib"
		fi
	],
	[
		lapack_lib=lapack
#		if   test "x$blas_lib" = "xatlas" ; then
#			lapack_lib=netlib
#		else
#			lapack_lib="$blaslib"
#		fi
	])

AC_ARG_WITH(mkl, [  --with-mkl=<path>       prefix of mkl installation. e.g. /op/intel/mkl81], [MKLDIR=$with_mkl], [MKLDIR=""])


AC_ARG_ENABLE([debug], [  --enable-debug          Turn on debugging],
	[
		case "${enableval}" in
		yes) debug=true ;;
		no)  debug=false ;;
		*) AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
		esac
	], [debug=false])

AC_ARG_ENABLE([gnuplot], [  --enable-gnuplot        Make use of gnuplot through pipes],
	[
		case "${enableval}" in
		yes) gnuplot=true ;;
		no)  gnuplot=false ;;
		*) AC_MSG_ERROR([bad value ${enableval} for --enable-gnuplot]) ;;
		esac
	], [gnuplot=true])

AC_ARG_ENABLE([mpi], [  --enable-mpi            Make use of mpi library],
	[
		case "${enableval}" in
		yes) use_mpi=true ;;
		no)  use_mpi=false ;;
		*) AC_MSG_ERROR([bad value ${enableval} for --enable-mpi]) ;;
		esac
	], [use_mpi=false])


# Check whether we are using pgplot
AC_ARG_ENABLE(pgplot, [  --enable-pgplot         Enable the pgplot library (default=no)],
	[
		if   test "x$enableval" = "xyes" ; then
			use_pgplot=yes
		else
			use_pgplot=no
		fi
	 ], [
		use_pgplot=no
	 ])

AM_CONDITIONAL([DEBUG], [test x$debug = xtrue])

if test "x$MKLDIR" != "" ; then
	AC_SUBST(MKLDIR)
	LDFLAGS="-L${MKLDIR}/lib/lib64 -L${MKLDIR}/lib/lib32 ${LDFLAGS}"
	CPPFLAGS="-I${MKLDIR}/include ${CPPFLAGS}"
fi


dnl -----------------------------------------------
dnl Check for libraries
dnl -----------------------------------------------


MPI_CXXLIBS=""
MPI_CXXFLAGS=""
use_mpi=no
#if test "x$use_mpi" = "xtrue" ; then
#	ACX_MPI([
#			use_mpi=yes
#			MPI_CXXLIBS=`mpicxx --showme:link`
#			MPI_CXXFLAGS=`mpicxx --showme:compile`
#			AC_DEFINE([HAVE_MPI], 1,[Define to 1 if you have mpi.])
#			AC_SUBST(HAVE_MPI)
#		],[
#			AC_MSG_ERROR([No mpi found. Configure without mpi, or adjust path])
#			use_mpi=no
#		])
#fi
AC_SUBST(MPI_CXXLIBS)
AC_SUBST(MPI_CXXFLAGS) 



AC_CHECK_LIB([m], [cos], [], AC_MSG_ERROR([Math library is required.]))

if test "x$blas_lib" = "xatlas" ; then
	AC_CHECK_LIB([atlas], [ATL_xerbla],[],[AC_CHECK_LIB([atlas_$ARCH],
		[ATL_xerbla], [], [AC_MSG_ERROR([Atlas library required!])])])
	AC_SEARCH_LIBS([cblas_dgemm], [atlas_$ARCH atlas cblas_$ARCH cblas], [], 
		[AC_MSG_ERROR([Atlas blas library is required.])])
	AC_DEFINE([HAVE_LIBCBLAS_ATLAS], 1,[Define to 1 if you have an ATLAS libcblas.])
	AC_SUBST(HAVE_LIBCBLAS_ATLAS)
	AC_DEFINE([HAVE_LIBCBLAS], 1,[Define to 1 if you have libcblas.])
	AC_SUBST(HAVE_LIBCBLAS)
elif test "x$blas_lib" = "xgsl" ; then
	AC_CHECK_LIB([gslcblas], [cblas_dgemm], [], AC_MSG_ERROR([GSL blas library required!]))
	AC_DEFINE([HAVE_LIBCBLAS_GSL], 1,[Define to 1 if you have an GSL libcblas.])
	AC_SUBST(HAVE_LIBCBLAS_GSL)
	AC_DEFINE([HAVE_LIBCBLAS], 1,[Define to 1 if you have libcblas.])
	AC_SUBST(HAVE_LIBCBLAS)
elif test "x$blas_lib" = "xcblas" ; then
	AC_CHECK_LIB([cblas], [cblas_dgemm], [], AC_MSG_ERROR([cblas library required!]))
	AC_DEFINE([HAVE_LIBCBLAS_CBLAS], 1,[Define to 1 if you have an GSL libcblas.])
	AC_SUBST(HAVE_LIBCBLAS_CBLAS)
	AC_DEFINE([HAVE_LIBCBLAS], 1,[Define to 1 if you have libcblas.])
	AC_SUBST(HAVE_LIBCBLAS)
elif test "x$blas_lib" = "xveclib" ; then
	AC_DEFINE([HAVE_LIBCBLAS_VECLIB], 1,[Define to 1 if you have a veclib framework.])
	AC_SUBST(HAVE_LIBCBLAS_VECLIB)
	AC_DEFINE([HAVE_LIBCBLAS], 1,[Define to 1 if you have libcblas.])
	AC_SUBST(HAVE_LIBCBLAS)
	LDFLAGS="-framework vecLib ${LDFLAGS}"
elif test "x$blas_lib" = "xmkl" ; then
	AC_CHECK_LIB([pthread], [pthread_create], [], AC_MSG_ERROR([GNU pthread library required!]))
	AC_CHECK_LIB([guide], [_end], [], AC_MSG_ERROR([Intel libguide library required!]))
	AC_SEARCH_LIBS([cblas_dgemm], [mkl_em64t mkl], [], [AC_MSG_ERROR([MKL blas library required!])])
	AC_DEFINE([HAVE_LIBCBLAS_MKL], 1,[Define to 1 if you have an MKL libcblas.])
	AC_SUBST(HAVE_LIBCBLAS_MKL)
	AC_DEFINE([HAVE_LIBCBLAS], 1,[Define to 1 if you have libcblas.])
	AC_SUBST(HAVE_LIBCBLAS)
fi


if test "x$lapack_lib" = "xnetlib" ; then
	AC_CHECK_HEADERS(f2c.h blaswrap.h)
	AC_SEARCH_LIBS([dsyevd_], [lapack_$ARCH lapack],
		[
		],[AC_MSG_ERROR([Netlib lapack library required!])])
	AC_SEARCH_LIBS([f2c_dgemm], [cblaswr_$ARCH cblaswr],[],
		[AC_MSG_ERROR([f2c (Fortran-to-C) library required!])])
	AC_SEARCH_LIBS([s_cat], [F77_$ARCH F77],[],
		[AC_MSG_ERROR([F77 (Fortran-77) library required!])])
	AC_DEFINE([HAVE_LIBLAPACK_NETLIB], 1,[Define to 1 if you have the netlib liblapack.])
	AC_SUBST(HAVE_LIBLAPACK_NETLIB)
	AC_DEFINE([HAVE_LIBLAPACK], 1,[Define to 1 if you have liblapack.])
	AC_SUBST(HAVE_LIBLAPACK)
elif test "x$lapack_lib" = "xclapack" ; then
	AC_CHECK_LIB([lapack], [dsyevd_], [], [AC_MSG_ERROR([The clapack library required!])])
	AC_DEFINE([HAVE_LIBLAPACK_CLAPACK], 1,[Define to 1 if you have the libclapack lapack library.])
	AC_SUBST(HAVE_LIBLAPACK_CLAPACK)
	AC_DEFINE([HAVE_LIBLAPACK], 1,[Define to 1 if you have liblapack.])
	AC_SUBST(HAVE_LIBLAPACK)
elif test "x$lapack_lib" = "xmkl" ; then
	if test "x$blas_lib" != "xmkl" ; then
		AC_CHECK_LIB([pthread], [pthread_create], [], AC_MSG_ERROR([GNU pthread library required!]))
		AC_CHECK_LIB([guide], [_end], [], AC_MSG_ERROR([Intel libguide library required!]))
	fi
	AC_CHECK_LIB([mkl_lapack], [dsyevd_], [],
		[AC_MSG_ERROR([MKL lapack library required!])])
	AC_DEFINE([HAVE_LIBLAPACK_MKL], 1,[Define to 1 if you have the libmkl lapack library.])
	AC_SUBST(HAVE_LIBLAPACK_MKL)
	AC_DEFINE([HAVE_LIBLAPACK], 1,[Define to 1 if you have a liblapack.])
	AC_SUBST(HAVE_LIBLAPACK)
elif test "x$lapack_lib" = "xgsl" ; then
	AC_DEFINE([HAVE_LIBLAPACK_GSL], 1,[Define to 1 if you only have the libgsl lapack library.])
	AC_SUBST(HAVE_LIBLAPACK_GSL)
elif test "x$lapack_lib" = "xlapack" ; then
	AC_CHECK_LIB([lapack], [dsyevd_], [],
		[AC_MSG_ERROR([Standard lapack library required!])])
	AC_DEFINE([HAVE_LIBLAPACK_LAPACK], 1,[Define to 1 if you have a standard lapack library.])
	AC_SUBST(HAVE_LIBLAPACK_LAPACK)
	AC_DEFINE([HAVE_LIBLAPACK], 1,[Define to 1 if you have a liblapack.])
	AC_SUBST(HAVE_LIBLAPACK)
elif test "x$lapack_lib" = "xveclib" ; then
	AC_DEFINE([HAVE_LIBLAPACK_VECLIB], 1,[Define to 1 if you have a veclib framework.])
	AC_SUBST(HAVE_LIBLAPACK_VECLIB)
	AC_DEFINE([HAVE_LIBLAPACK], 1,[Define to 1 if you have a liblapack.])
	AC_SUBST(HAVE_LIBLAPACK)
fi


# We also need pgplot for the tempo2 plugin
AM_CONDITIONAL([PGPLOT], [test x$use_pgplot = xyes])
AM_CONDITIONAL([PGPLOT2], [test x$use_pgplot = xyes])
if test "x$use_pgplot" != "xno" ; then
	AC_CHECK_LIB([z],       [deflate],   [],
		[], [])
	AC_CHECK_HEADERS(cpgplot.h)
	AC_CHECK_LIB([png],     [png_error], [],
		[AC_MSG_ERROR([png library is required])], [])
	AC_CHECK_LIB(X11, XOpenDisplay, , [
			echo "ERROR:  You need libX11 to build tempo2.  Verify that you have libX11.a or";
			echo "        libX11.so installed and that it is located in the X libraries";
			echo "        directory shown above.  If it is in a different directory, try using";
			echo "        the --x-libraries parameter to configure.";
		AC_MSG_ERROR([Fatal:  libX11 not found.])
		])
	AC_CHECK_LIB([g2c], [do_fio],,
		[
			AC_CHECK_LIB([gfortran], [open],,
				[
					AC_MSG_ERROR([Cannot link to libg2c or libfortran])
				])
		])

	AC_CHECK_LIB([pgplot],  [pgbeg_],    [],
		[
			AC_MSG_ERROR([pgplot library is required])
		], [])
	AC_CHECK_LIB([cpgplot], [cpgbeg],    [],
		[AC_MSG_ERROR([pgplot library is required])], [])

	AC_DEFINE([HAVE_PGPLOT], 1,[Define to 1 if you have the libpgplot library.])
	AC_SUBST(HAVE_PGPLOT)
fi

if test "x$gnuplot" != "xno" ; then
	AC_DEFINE([HAVE_GNUPLOT], 1,[Define to 1 if you have gnuplot installed.])
	AC_SUBST(HAVE_GNUPLOT)
fi


AC_CHECK_LIB([gsl], [gsl_acosh], [], AC_MSG_ERROR([GNU Scientific library is required.]))


# Registering the variables so they are known in the Makefile's
AC_SUBST(DEPS_CFLAGS)
AC_SUBST(DEPS_LIBS)

AC_CONFIG_MACRO_DIR([m4])

dnl -----------------------------------------------
dnl Generates Makefile's, configuration files and scripts
dnl -----------------------------------------------

# Show the scripts where the Makefile's are
AC_CONFIG_FILES([Makefile src/Makefile src/tools/Makefile])
AC_CONFIG_COMMANDS([default],[[ echo timestamp > stamp-h ]],[[]])
AC_OUTPUT
