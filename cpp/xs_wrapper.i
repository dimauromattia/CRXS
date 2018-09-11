%module xs_wrapper
%include "std_string.i"

%{
    #define SWIG_FILE_WITH_INIT
    #include "xs_wrapper.h"
%}


%include "xs_wrapper.h"

