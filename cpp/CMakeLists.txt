cmake_minimum_required(VERSION 2.8)
project(CRXS)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-std=c++11)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)

set(BIN     ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(LIB     ${CMAKE_CURRENT_SOURCE_DIR}/lib)

message(STATUS Set\ bin\ directory\ to:\  ${BIN}    )
message(STATUS Set\ lib\ directory\ to:\  ${LIB}    )
message(STATUS Copy\ include\ files\ to:\ ${INCLUDE})

find_package(GSL REQUIRED)

include_directories(${INCLUDE})
include_directories(${GSL_INCLUDE_DIRS}/gsl)

message(STATUS GSL\ include\ files\ to:\ ${GSL_INCLUDE_DIRS}/gsl)



#OwnLibraries
link_directories(${LIB})
IF(CMAKE_SYSTEM_NAME STREQUAL Darwin)
  set(myLIB ${CMAKE_CURRENT_SOURCE_DIR}/lib/libCRXS.dylib)
ELSEIF(CMAKE_SYSTEM_NAME STREQUAL Linux)
  set(myLIB ${CMAKE_CURRENT_SOURCE_DIR}/lib/libCRXS.so)
ENDIF()


add_subdirectory(src)
