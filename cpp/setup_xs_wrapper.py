#!/usr/bin/env python3 

"""
setup.py For pbar XS software
"""

from distutils.core import setup, Extension
import numpy
import os

try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()
    
CDIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)) )

example_module = Extension('_xs_wrapper',
                           sources=['xs_wrapper_wrap.cxx', 'xs_wrapper.cxx'],
                           include_dirs=[CDIR+'/cpp/include',numpy_include],
                           library_dirs=[CDIR+'/cpp/lib'],
                           libraries=['CRXS'],
                           extra_compile_args=['-std=c++11']
                           )


setup (name = 'xs_wrapper',
       version = '1.0',
       author      = "Michael Korsmeier",
       description = """CR XS software""",
       ext_modules = [example_module],
       py_modules = ["xs_wrapper"],
       )
