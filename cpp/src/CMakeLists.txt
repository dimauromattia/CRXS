#Library
add_library(CRXS SHARED xs_definitions.cxx
                        xs_definitions.h
                        xs.cxx
                        xs_pbar.cxx
                        xs_Dbar.cxx
                        xs.h                )



target_link_libraries(CRXS GSL::gsl GSL::gslcblas)

file(  COPY xs_definitions.h DESTINATION ${INCLUDE}  )
file(  COPY xs.h             DESTINATION ${INCLUDE}  )


add_subdirectory(example)
