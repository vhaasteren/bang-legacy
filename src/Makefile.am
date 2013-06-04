bindir = $(exec_prefix)/bin
libdir = $(exec_prefix)/lib

ban_sources = corefunctions.cpp analyze.cpp linal.cpp banfunc.cpp linalfunc.cpp filefunctions.cpp correlation.cpp moremath.cpp
ban_headers = corefunctions.h analyze.h linal.h banfunc.h linalfunc.h filefunctions.h moremath.h correlation.h

ban_gui_sources = bangui.cpp
ban_gui_headers = bangui.h

iugui_sources = iugui.cpp
iugui_headers = iugui.h

ban_libs = libbanginterface.la

PTALIBS = $(ban_libs)
PTAHEADERS = $(ban_headers) $(ban_gui_headers) $(iugui_headers)

# These are the main program
bin_PROGRAMS =
lib_LTLIBRARIES = $(PTALIBS)
noinst_HEADERS = $(PTAHEADERS)

# This is the interface exposed to python/C
libbanginterface_la_SOURCES = banginterface.cpp $(ban_gui_sources) $(iugui_sources) $(ban_sources)
libbanginterface_la_LIBADD = $(DEPS_LIBS)
libbanginterface_la_LDFLAGS = -version-info $(GENERIC_LIBRARY_VERSION) -release $(GENERIC_RELEASE)

AM_CPPFLAGS = $(DEPS_CFLAGS) $(MPI_CXXFLAGS)

SUBDIRS = tools 
