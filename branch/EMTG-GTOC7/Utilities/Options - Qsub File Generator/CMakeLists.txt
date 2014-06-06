cmake_minimum_required(VERSION 2.4)


option(TAUB_BOOST_LOCATION "Set that you want to use the TAUB Boost Location
Specifically, otherwise, assume defaults work." OFF)



if(TAUB_BOOST_LOCATION)
message(STATUS "We are on Taub.  Forcing to /usr/local/boost/boost-1.51.0")
	set(BOOST_ROOT "/usr/local/boost/boost-1.51.0")
	set(Boost_NO_SYSTEM_PATHS ON)
else(TAUB_BOOST_LOCATION)
	set(BOOST_ROOT ../../boost) #set a hint to find boost, but it'll search in the obvious spots too
	set(BOOST_INCLUDEDIR ${BOOST_ROOT}/boost)
	set(BOOST_LIBRARYDIR ${BOOST_ROOT}/stage/lib)
endif(TAUB_BOOST_LOCATION)

project(QsubFileGenerator)
include(FindBoost)
find_package(Boost 1.44 COMPONENTS filesystem system REQUIRED) #minimum version is 1.44, REQUIRED means we fail cmake if we can't find boost

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

add_subdirectory(${PROJECT_SOURCE_DIR}/src)
