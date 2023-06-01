%module xs_wrapper
%include "std_string.i"

%{
    #define SWIG_FILE_WITH_INIT
    #include "xs_wrapper.h"
    #include "string"
%}

%include "numpy.i"

%init %{
    import_array();
%}

%apply (double* IN_ARRAY1, int DIM1) {(double* C_array, int len_C_array)};
%apply (double* IN_ARRAY1, int DIM1) {(double* D_array, int len_D_array)};

%include "xs_wrapper.h"

