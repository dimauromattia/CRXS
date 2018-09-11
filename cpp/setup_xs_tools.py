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


module         = Extension('_xs_tools',
                           sources=['xs_tools_wrap.cxx', 'xs_tools.cxx'],
                           include_dirs=[CDIR+'/cpp/include', numpy_include],
                           library_dirs=[CDIR+'/cpp/lib'],
                           libraries=['CRXS'],
                           extra_compile_args=['-std=c++11']
                           )

setup (name = 'xs_tools',
       version = '1.0',
       author      = "Michael Korsmeier",
       description = """CRXS software""",
       ext_modules = [module],
       py_modules = ["xs_tools"],
       )
