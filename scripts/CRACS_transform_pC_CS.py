#! /usr/bin/env python

import os
import sys
import glob
import math
import numpy as np


CRACS = os.environ['CRACS']
data_dir    = CRACS+'/data/CS_measurements/pC'
res_dir     = data_dir
script_dir  = CRACS+'/scripts/CS_transformation_pC'

sys.path.append( script_dir )


import na49
na49.write( data_dir+'/na49',  res_dir )

import abramov
abramov.write( data_dir+'/abramov',  res_dir )

import amann
amann.write( data_dir+'/amann',  res_dir )

import barton
barton.write( data_dir+'/barton',  res_dir )

import sugaya
sugaya.write( data_dir+'/sugaya',  res_dir )
