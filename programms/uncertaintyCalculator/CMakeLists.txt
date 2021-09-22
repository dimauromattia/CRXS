add_executable(CRACS_uncertaintyCalculator uncertaintyCalculator.cpp)

target_link_libraries(CRACS_uncertaintyCalculator  Qt5::Gui Qt5::Xml)
target_link_libraries(CRACS_uncertaintyCalculator  ${myLIB})
target_link_libraries(CRACS_uncertaintyCalculator  ${ROOTLIBS})
#target_link_libraries(CRACS_uncertaintyCalculator  ${MULTINESTLIB_MPI})
#target_link_libraries(CRACS_uncertaintyCalculator  ${MPI_CXX_LIBRARIES})

add_executable(CRACS_uncertaintyCalculator_old uncertaintyCalculator_old.cpp)

target_link_libraries(CRACS_uncertaintyCalculator_old  Qt5::Gui Qt5::Xml)
target_link_libraries(CRACS_uncertaintyCalculator_old  ${myLIB})
target_link_libraries(CRACS_uncertaintyCalculator_old  ${ROOTLIBS})
#target_link_libraries(CRACS_uncertaintyCalculator_old  ${MULTINESTLIB_MPI})
#target_link_libraries(CRACS_uncertaintyCalculator_old  ${MPI_CXX_LIBRARIES})

