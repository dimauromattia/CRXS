add_executable(CRACS_crossSectionFitter crossSectionFitter.cpp)

target_link_libraries(CRACS_crossSectionFitter  Qt5::Gui Qt5::Xml)
target_link_libraries(CRACS_crossSectionFitter  ${myLIB})
target_link_libraries(CRACS_crossSectionFitter  ${ROOTLIBS})


