add_executable(CRXS_test test.cpp)
target_link_libraries(CRXS_test ${myLIB})
target_link_libraries(CRXS_test GSL::gsl GSL::gslcblas)


