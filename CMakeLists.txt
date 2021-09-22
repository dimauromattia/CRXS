cmake_minimum_required(VERSION 2.8)
project(CRACSproject)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")
message("CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)


set(BIN ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(INCLUDEDIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(LIB ${CMAKE_CURRENT_SOURCE_DIR}/lib)

message(STATUS Set\ bin\ directory\ to:\ ${BIN})
message(STATUS Set\ lib\ directory\ to:\ ${LIB})
message(STATUS Copy\ include\ files\ to:\ ${INCLUDEDIR})


include_directories(${INCLUDEDIR})

#Add root
include_directories($ENV{ROOTSYS}/include			)
link_directories   ($ENV{ROOTSYS}/lib				)
set(ROOTLIBS 		$ENV{ROOTSYS}/lib/libCore.so 
					$ENV{ROOTSYS}/lib/libHist.so
					$ENV{ROOTSYS}/lib/libMathCore.so
 					$ENV{ROOTSYS}/lib/libMatrix.so
 					$ENV{ROOTSYS}/lib/libRIO.so
 					$ENV{ROOTSYS}/lib/libThread.so
  					$ENV{ROOTSYS}/lib/libGraf.so
 					$ENV{ROOTSYS}/lib/libTree.so
 					$ENV{ROOTSYS}/lib/libGpad.so
 					$ENV{ROOTSYS}/lib/libMinuit.so
 					$ENV{ROOTSYS}/lib/libMinuit2.so
 					$ENV{ROOTSYS}/lib/libGraf3d.so   )
message(STATUS ${ROOTLIBS})

#Add QT
set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(Qt5Gui REQUIRED)
find_package(Qt5Xml REQUIRED)

#Add cfitsio
#set(CFITSIOLIB 			$ENV{CFITSIO_DIR}/lib/libcfitsio.a	)
#include_directories   (	$ENV{CFITSIO_DIR}/include			)
#link_directories      (	$ENV{CFITSIO_DIR}/lib         		)

#Add multinest
IF(CMAKE_SYSTEM_NAME STREQUAL Darwin)
	set(MULTINESTLIB 		$ENV{MULTINEST_DIR}/lib/libmultinest.dylib		)
	set(MULTINESTLIB_MPI 	$ENV{MULTINEST_DIR}/lib/libmultinest_mpi.dylib	)
ELSEIF(CMAKE_SYSTEM_NAME STREQUAL Linux)
	set(MULTINESTLIB		$ENV{MULTINEST_DIR}/lib/libmultinest.so			)
	set(MULTINESTLIB_MPI 	$ENV{MULTINEST_DIR}/lib/libmultinest_mpi.so		)
ENDIF()
include_directories		(   $ENV{MULTINEST_DIR}/include )
link_directories		(   $ENV{MULTINEST_DIR}/lib		)

#Add mpi

#find_package(MPI)
#message(STATUS Use\ the\ follwoing\ mpi\ compile\ flags:\  ${MPI_CXX_LIBRARIES})
#message(STATUS Use\ the\ follwoing\ mpi\ compile\ flags:\  ${MPI_INCLUDE_PATH})

#include_directories(SYSTEM ${MPI_INCLUDE_PATH})

#OwnLibraries
link_directories(${LIB})
IF(CMAKE_SYSTEM_NAME STREQUAL Darwin)
    set(myLIB ${CMAKE_CURRENT_SOURCE_DIR}/lib/libCRACS.dylib)
ELSEIF(CMAKE_SYSTEM_NAME STREQUAL Linux)
    set(myLIB ${CMAKE_CURRENT_SOURCE_DIR}/lib/libCRACS.so)
ENDIF()


#Add libraries
add_subdirectory(libraries)

#Add programms
add_subdirectory(programms)



