cmake_minimum_required(VERSION 2.4)


option(TAUB_BOOST_LOCATION "Set that you want to use the TAUB Boost Location
Specifically, otherwise, assume defaults work." OFF)
option(BLUE_WATERS "Is this on blue waters?" OFF)

if(BLUE_WATERS)
	message(STATUS "We are on Blue Waters.")
	set(CMAKE_CXX_COMPILER CC)
	set(CMAKE_C_COMPILER cc)
	execute_process(COMMAND "module unload gcc" COMMAND "module load gcc/4.7.1")
	set(BOOST_ROOT "/sw/xe/boost/1.53.0/sles11.1_gnu4.8.1") 
elseif(TAUB_BOOST_LOCATION)
	message(STATUS "We are on Taub.  Forcing to /usr/local/boost/boost-1.51.0")
	set(BOOST_ROOT "/usr/local/boost/boost-1.51.0")
	set(Boost_NO_SYSTEM_PATHS ON)
elseif(MSVC)
	set(BOOST_ROOT ../../boost) #set a hint to find boost, but it'll search in the obvious spots too
	set(BOOST_INCLUDEDIR ${BOOST_ROOT}/boost)
	set(BOOST_LIBRARYDIR ${BOOST_ROOT}/stage/lib)
endif(BLUE_WATERS)


project(QsubFileGenerator)

include(FindBoost)
find_package(Boost 1.44 COMPONENTS filesystem system REQUIRED) #minimum version is 1.44, REQUIRED means we fail cmake if we can't find boost

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

add_subdirectory(${PROJECT_SOURCE_DIR}/src)
