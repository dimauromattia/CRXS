add_executable(CRACS_lisFluxFitter lisFluxFitter.cpp)

target_link_libraries(CRACS_lisFluxFitter  Qt5::Gui Qt5::Xml)
target_link_libraries(CRACS_lisFluxFitter  ${myLIB})
target_link_libraries(CRACS_lisFluxFitter  ${ROOTLIBS})
#target_link_libraries(CRACS_lisFluxFitter  ${MULTINESTLIB_MPI})
#target_link_libraries(CRACS_lisFluxFitter  ${MPI_CXX_LIBRARIES})

