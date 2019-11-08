#Library

configure_file(crxs.cpp.in crxs.cpp @ONLY)

add_library(CRXS SHARED crxs.h
                        ${CMAKE_CURRENT_BINARY_DIR}/crxs.cpp
                        xs_definitions.cxx
                        xs_definitions.h
                        xs.cxx
                        xs_pbar.cxx
                        xs_p.cxx
                        xs_Dbar.cxx
                        xs.h
                        linAlg_tools.cxx
                        linAlg_tools.h           )


target_link_libraries(CRXS GSL::gsl GSL::gslcblas)

file(  COPY xs_definitions.h    DESTINATION ${INCLUDE}  )
file(  COPY xs.h                DESTINATION ${INCLUDE}  )
file(  COPY crxs.h              DESTINATION ${INCLUDE}  )
file(  COPY linAlg_tools.h      DESTINATION ${INCLUDE}  )

#add_subdirectory(example)
