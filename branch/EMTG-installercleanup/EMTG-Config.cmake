#EMTG Cmake pre-build configuration file.
#This file contains some pre-defined addresses that you may need to customize for your own system.  Setting the 
#below correctly can assist with building EMTG quickly and easily.


#---------SNOPT hints------------
#Change this next line to point towards your snopt7 root directory.  This is not your cppsrc or your lib directory, but up one level from that.
#If SNOPT has been installed on the system path, this hint is probably unnecessary. 
#You can alternatively add and set the SNOPT_INCLUDEDIR and SNOPT_LIBDIR as direct paths to your cppsrc and appropriate lib directory.

#uncomment this line and change it to be your appropriate path (works for Windows or Unix-based systems)
	#set(SNOPT_ROOT_DIR /home/user/snopt)
	#set(SNOPT_INCLUDEDIR /home/user/snopt/cppsrc)
	#set(SNOPT_LIBDIR /home/user/snopt/lib)

#-------BOOST HINTS------------------
#cmake usually has an easy time finding a properly installed boost.  If it cannot find your copy of boost, uncomment the next two lines and 
#appropriately modify them to point at your boost distribution.

	#set(BOOST_ROOT ../../boost) #set a hint to find boost, but it'll search in the obvious spots too
	#set(BOOST_INCLUDEDIR ${BOOST_ROOT}/boost) #frequently this is where the headers are
	#set(BOOST_LIBRARYDIR ${BOOST_ROOT}/lib)   # this is sometimes {BOOST_ROOT}/stage/lib	
	
	#This option will tell boost NOT to look at the system path.  If you want to force boost to use some local version, this is how
	#set(Boost_NO_SYSTEM_PATHS ON)
	
	
#----------MPI SNOPT------------------
#If you are going to build SNOPT to operate on a system using MPI, this must be set in this config file before
#anything else gets run.  Uncomment the next three lines, and change the last two appropriately.

	#set(EMTG_MPI ON)
	#set(CMAKE_CXX_COMPILER mpicxx) #set 'mpicxx' to your MPI distribution's chosen c++ compiler alias that is on your path
	#set(CMAKE_C_COMPILER mpicc) #set 'mpicc' to be your MPI distribution's chosen c++ compiler alias that is on your path.

#Note to windows users: if you've built and MPI server in windows and properly set it up to build within and against Visual studio projects
#the bottom two lines above are NOT required, but the EMTG_MPI line still must be uncommented.
