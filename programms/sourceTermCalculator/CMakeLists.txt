add_executable(CRACS_sourceTermCalculator sourceTermCalculator.cpp)

target_link_libraries(CRACS_sourceTermCalculator  Qt5::Gui Qt5::Xml)
target_link_libraries(CRACS_sourceTermCalculator  ${myLIB})
target_link_libraries(CRACS_sourceTermCalculator  ${ROOTLIBS})
#target_link_libraries(CRACS_sourceTermCalculator  ${MULTINESTLIB_MPI})
#target_link_libraries(CRACS_sourceTermCalculator  ${MPI_CXX_LIBRARIES})

