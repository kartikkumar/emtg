cmake_minimum_required(VERSION 2.4)

#set(CMAKE_C_COMPILER mpicc)
#set(CMAKE_CXX_COMPILER mpic++)

#identify the compiler
project(EMTG)



#-------------do we support C++11?--------------------------------
message(STATUS "Checking C++11 support...")
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-std=c++11")
	message(STATUS "Full support of C++11 found.")
elseif(COMPILER_SUPPORTS_CXX0X)
	set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-std=c++0x")
	message(STATUS "Support for C++0x found.")
else(COMPILER_SUPPORTS_CXX11)
	add_definitions(-D_STONEAGECplusplus)
	message(STATUS "Stoneage C++ mode enabled.")
endif(COMPILER_SUPPORTS_CXX11)




#compiler-time Define-level options  for EMTG
#option(EMTG_MPI "Enable the parallel functions of EMTG" OFF) // this one pretty much has to be at cmake run-time due to the compiler differences
option(QUIET_SNOPT "If you have a special build of SNOPT, this will silence it completely." OFF)
option(BACKGROUND_MODE "This suppresses 'press any key to continue' at the end of an EMTG run" OFF)
option(EMTG_proprietary_access "This turns on Goddard-internal proprietary hooks for code not distributed with standard EMTG" OFF)

mark_as_advanced(QUIET_SNOPT EMTG_proprietary_access) #these two are only useful if you are fully aware of what they do, hence hide them from cmake-gui/ccmake

if(EMTG_MPI)
	add_definitions(-DEMTG_MPI)
endif(EMTG_MPI)

if(QUIET_SNOPT)
	add_definitions(-DQUIET_SNOPT)
endif(QUIET_SNOPT)

if(BACKGROUND_MODE)
	add_definitions(-DBACKGROUND_MODE)
endif(BACKGROUND_MODE)

if(EMTG_proprietary_access)
	add_definitions(-D_EMTG_proprietary)
endif(EMTG_proprietary_access)





#----------------CSpice Work------------------------------------------------------
message(STATUS "Now checking for CSpice...")
if (MSVC) # aka 'windows'
	set(CSPICE_FILE_TO_BE_CORRECTED "${CMAKE_SOURCE_DIR}/cspice/src/cspice/mkprodct.bat")
	set(CSPICE_MAKEFILE "makeall.bat")
else (MSVC) #aka Linux or Mac, as both of those use the mkprodct.csh system
	set(CSPICE_FILE_TO_BE_CORRECTED "${CMAKE_SOURCE_DIR}/cspice/src/cspice/mkprodct.csh")
	set(CSPICE_MAKEFILE "makeall.csh")
endif(MSVC)


if (EXISTS "${CMAKE_SOURCE_DIR}/cspice/${CSPICE_MAKEFILE}")
	message(STATUS "Cspice located and unzipped correctly.")
else()
	message(FATAL_ERROR "CSPICE is not unzipped into the emtg/cspice folder correctly.")
endif()

if (NOT EXISTS "${CMAKE_SOURCE_DIR}/cspice/lib/libcspice.a" AND NOT EXISTS "${CMAKE_SOURCE_DIR}/cspice/lib/cspice.lib") # we need to build cspice

	message(STATUS "CSpice needs to be built.... adding -DMIX_C_AND_FORTRAN flag and attempting build.")	


	if (MSVC)
		file(GLOB_RECURSE FILE_LIST "${CMAKE_SOURCE_DIR}/cspice/src/*.bat")
		foreach(CSPICE_FILE_TO_BE_CORRECTED ${FILE_LIST})
			FILE(READ ${CSPICE_FILE_TO_BE_CORRECTED} FILE_CONTENT)
			STRING(REPLACE "-DNON_ANSI_STDIO" "-DNON_ANSI_STDIO -DMIX_C_AND_FORTRAN" MODIFIED_FILE_CONTENT "${FILE_CONTENT}")
			FILE(WRITE ${CSPICE_FILE_TO_BE_CORRECTED} "${MODIFIED_FILE_CONTENT}")
		endforeach()
	else(MSVC)
		file(GLOB_RECURSE FILE_LIST "${CMAKE_SOURCE_DIR}/cspice/src/*.csh")
		 foreach(CSPICE_FILE_TO_BE_CORRECTED ${FILE_LIST})
                        FILE(READ ${CSPICE_FILE_TO_BE_CORRECTED} FILE_CONTENT)
			STRING(REPLACE "-DNON_UNIX_STDIO" "-DNON_UNIX_STDIO -DMIX_C_AND_FORTRAN" MODIFIED_FILE_CONTENT "${FILE_CONTENT}")                    
                        FILE(WRITE ${CSPICE_FILE_TO_BE_CORRECTED} "${MODIFIED_FILE_CONTENT}")
                endforeach()

	endif(MSVC)

	execute_process(COMMAND "${CMAKE_SOURCE_DIR}/cspice/${CSPICE_MAKEFILE}" WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/cspice")

	if (NOT MSVC) #cspice doesn't label the file correctly when it is done; easier to just rename it
		message(STATUS "Renaming cspice.a to libcspice.a")
		file(RENAME "${CMAKE_SOURCE_DIR}/cspice/lib/cspice.a" "${CMAKE_SOURCE_DIR}/cspice/lib/libcspice.a")
		if (NOT EXISTS "${CMAKE_SOURCE_DIR}/cspice/lib/libcspice.a")
			message(FATAL_ERROR "CSpice not built properly for a Unix System.")
		endif()
	else()
		 if (NOT EXISTS "${CMAKE_SOURCE_DIR}/cspice/lib/cspice.lib")
                        message(FATAL_ERROR "CSpice not built properly for a Unix System.")
                endif()
	endif(NOT MSVC)
endif()

set(EMTG_LIBRARIES ${EMTG_LIBRARIES} cspice)
include_directories("${CMAKE_SOURCE_DIR}/cspice/include")
link_directories("${CMAKE_SOURCE_DIR}/cspice/lib")

#----------------------------SNOPT Stuff-------------------------------------
message(STATUS "Now checking on SNOPT...")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR})
include(FindSNOPT)
find_package(SNOPT REQUIRED)

include_directories(${SNOPT_INCLUDE_DIR})
set(EMTG_LIBRARIES ${EMTG_LIBRARIES} ${SNOPT_LIBRARIES})


#-------------------------Boost related work------------------------------
message(STATUS "Now looking for BOOST..")
include(FindBoost)
find_package(Boost 1.44 COMPONENTS filesystem system serialization REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

set(EMTG_LIBRARIES ${EMTG_LIBRARIES} ${Boost_LIBRARIES})

#pull in EMTG itself
add_subdirectory(src)
