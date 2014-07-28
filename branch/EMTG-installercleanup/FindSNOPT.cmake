# Find SNOPT
#
# This sets the following variables:
# SNOPT_FOUND
# SNOPT_INCLUDE_DIRS
# SNOPT_LIBRARIES
# SNOPT_DEFINITIONS

find_package(PkgConfig QUIET)
pkg_check_modules(SNOPT snopt)
set(SNOPT_DEFINITIONS ${SNOPT_CFLAGS_OTHER})

#this will find the directory with the snoptProblem.h.  It should also be the same directory we need for the rest of the headers.  Whether old SNOPT or new SNOPT, all headers are roughly the same and in the same directory
find_path(SNOPT_INCLUDE_DIR cppsrc/snoptProblem.hh cppsrc/snoptProblem.hpp snoptProblem.hh snoptProblem.hpp
    HINTS ${SNOPT_INCLUDEDIR} ${SNOPT_INCLUDE_DIRS}
    PATHS "${CMAKE_INSTALL_PREFIX}/include" "${PROJECT_SOURCE_DIR}/SNOPT/cppsrc" "${CMAKE_SOURCE_DIR}/SNOPT/interfaces/cppsrc")


#now we need to look for the snopt libs and the lib directory. Depending on if we are heritage or not, we need to find and collect different sets of libraries.  We'll use the abse snopt library as the indicator

find_library(SNOPT_PRIMARY_LIBRARY NAMES snopt_c
             HINTS ${SNOPT_LIBDIR} ${SNOPT_LIBRARY_DIRS}
	     PATHS "${PROJECT_SOURCE_DIR}/SNOPT/lib")

if(SNOPT_PRIMARY_LIBRARY-NOTFOUND) #must be the new library

	find_library(SNOPT_PRIMARY_LIBRARY NAMES snopt7
             HINTS ${SNOPT_LIBDIR} ${SNOPT_LIBRARY_DIRS}
	     PATHS "${PROJECT_SOURCE_DIR}/SNOPT/lib/.libs")

	find_library(SNOPT_INTERFACE NAMES snopt7_cpp  HINTS ${SNOPT_LIBDIR} ${SNOPT_LIBRARY_DIRS}
             PATHS "${PROJECT_SOURCE_DIR}/SNOPT/lib/.libs")

else(SNOPT_PRIMARY_LIBRARY-NOTFOUND) #must still be the old interface

 find_library(SNOPT_INTERFACE NAMES snopt_cpp  HINTS ${SNOPT_LIBDIR} ${SNOPT_LIBRARY_DIRS}
             PATHS "${PROJECT_SOURCE_DIR}/SNOPT/lib")

# find_library(SNOPT_PRINT NAMES snprint_c  HINTS ${SNOPT_LIBDIR} ${SNOPT_LIBRARY_DIRS} PATHS "${PROJECT_SOURCE_DIR}/SNOPT/lib")

 find_library(SNOPT_F2C NAMES f2c  HINTS ${SNOPT_LIBDIR} ${SNOPT_LIBRARY_DIRS}
             PATHS "${PROJECT_SOURCE_DIR}/SNOPT/lib" NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

	add_definitions(-DHeritage_SNOPT7)
endif(SNOPT_PRIMARY_LIBRARY-NOTFOUND)

set(SNOPT_LIBRARY ${SNOPT_LIBRARIES} ${SNOPT_PRIMARY_LIBRARY} ${SNOPT_INTERFACE} ${SNOPT_F2C})
set(SNOPT_LIBRARIES ${SNOPT_LIBRARY})

set(SNOPT_INCLUDE_DIRS ${SNOPT_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SNOPT DEFAULT_MSG
                                  SNOPT_LIBRARY SNOPT_INCLUDE_DIR)

mark_as_advanced(SNOPT_INCLUDE_DIRS SNOPT_LIBRARIES)
