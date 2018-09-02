#!/usr/bin/env python3 

"""
setup.py For pbar XS software
"""

from distutils.core import setup, Extension
import numpy
import os

ROOTSYS = os.environ['ROOTSYS']
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()


example_module = Extension('_clike',
                           sources=['clike_wrap.cxx', 'clike.cxx'],
                           include_dirs=[numpy_include],
                           )

setup (name = 'clike',
       version = '1.0',
       author      = "Michael Korsmeier",
       description = """pbar XS software""",
       ext_modules = [example_module],
       py_modules = ["clike"],
       )
