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
    PATHS "${CMAKE_INSTALL_PREFIX}/include" "${PROJECT_SOURCE_DIR}/SNOPT/cppsrc" "${CMAKE_SOURCE_DIR}/SNOPT/interfaces/cppsrc" "${SNOPT_INCLUDEDIR}" "${SNOPT_ROOT_DIR}/cppsrc" "${SNOPT_ROOT_DIR}/interfaces/cppsrc")


#now we need to look for the snopt libs and the lib directory. Depending on if we are heritage or not, we need to find and collect different sets of libraries.  We'll use the abse snopt library as the indicator

find_library(SNOPT_PRIMARY_LIBRARY NAMES snopt_c
             HINTS ${SNOPT_LIBDIR} ${SNOPT_LIBRARY_DIRS}
	     PATHS "${PROJECT_SOURCE_DIR}/SNOPT/lib" "${SNOPT_LIBDIR}" "${SNOPT_ROOT_DIR}/lib")
	if (NOT SNOPT_PRIMARY_LIBRARY)
		set(NEW_SNOPT 1 CACHE BOOL "save new SNOPT flag")
	endif(NOT SNOPT_PRIMARY_LIBRARY)
		 
if(NEW_SNOPT) #must be the new library
	find_library(SNOPT_PRIMARY_LIBRARY NAMES snopt7
             HINTS ${SNOPT_LIBDIR} ${SNOPT_LIBRARY_DIRS}
	     PATHS "${PROJECT_SOURCE_DIR}/SNOPT/lib/.libs" "${SNOPT_LIB}" "${SNOPT_ROOT_DIR}/lib/.libs" "${SNOPT_ROOT_DIR}/lib")

	find_library(SNOPT_INTERFACE NAMES snopt7_cpp  HINTS ${SNOPT_LIBDIR} ${SNOPT_LIBRARY_DIRS}
             PATHS "${PROJECT_SOURCE_DIR}/SNOPT/lib/.libs" "${SNOPT_LIB}" "${SNOPT_ROOT_DIR}/lib/.libs" "${SNOPT_ROOT_DIR}/lib")
	set(SNOPT_INCLUDE_DIR ${SNOPT_INCLUDE_DIR} ${SNOPT_INCLUDE_DIR}/..)
	set(SNOPT_INCLUDE_DIRS ${SNOPT_INCLUDE_DIR})
	
	set(SNOPT_LIBRARY ${SNOPT_LIBRARIES} ${SNOPT_PRIMARY_LIBRARY} ${SNOPT_INTERFACE})
else(NEW_SNOPT) #must still be the old interface

 find_library(SNOPT_INTERFACE NAMES snopt_cpp  HINTS ${SNOPT_LIBDIR} ${SNOPT_LIBRARY_DIRS}
             PATHS "${PROJECT_SOURCE_DIR}/SNOPT/lib"  "${SNOPT_LIBDIR}" "${SNOPT_ROOT_DIR}/lib")

 find_library(SNOPT_F2C NAMES f2c  HINTS ${SNOPT_LIBDIR} ${SNOPT_LIBRARY_DIRS}
             PATHS "${PROJECT_SOURCE_DIR}/SNOPT/lib" "${SNOPT_LIBDIR}" "${SNOPT_ROOT_DIR}/lib" NO_SYSTEM_ENVIRONMENT_PATH NO_CMAKE_SYSTEM_PATH)

	add_definitions(-DHeritage_SNOPT7)
	set(SNOPT_LIBRARY ${SNOPT_LIBRARIES} ${SNOPT_PRIMARY_LIBRARY} ${SNOPT_INTERFACE} ${SNOPT_F2C})
endif(NEW_SNOPT)


set(SNOPT_LIBRARIES ${SNOPT_LIBRARY})

set(SNOPT_INCLUDE_DIRS ${SNOPT_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SNOPT DEFAULT_MSG
                                  SNOPT_LIBRARY SNOPT_INCLUDE_DIR)

mark_as_advanced(SNOPT_INCLUDE_DIRS SNOPT_LIBRARIES)
