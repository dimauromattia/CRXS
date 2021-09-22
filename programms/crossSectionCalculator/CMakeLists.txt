add_executable(CRACS_crossSectionCalculator crossSectionCalculator.cpp)

target_link_libraries(CRACS_crossSectionCalculator  Qt5::Gui Qt5::Xml)
target_link_libraries(CRACS_crossSectionCalculator  ${myLIB})
#target_link_libraries(CRACS_crossSectionCalculator  ${ROOTLIBS})
#target_link_libraries(CRACS_crossSectionCalculator  ${MULTINESTLIB_MPI})
#target_link_libraries(CRACS_crossSectionCalculator  ${MPI_CXX_LIBRARIES})

