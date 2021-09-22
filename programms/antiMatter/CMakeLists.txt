add_executable(CRACS_AntiMatter antiMatter.cpp)

target_link_libraries(CRACS_AntiMatter  Qt5::Gui Qt5::Xml)
target_link_libraries(CRACS_AntiMatter  ${myLIB})
target_link_libraries(CRACS_AntiMatter  ${ROOTLIBS})
