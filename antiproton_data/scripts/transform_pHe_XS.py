#! /usr/bin/env python

import os
import sys
import glob
import math
import numpy as np



CDIR = os.path.dirname(os.path.realpath(__file__))

data_dir = CDIR+'/../data/pHe'
res_dir  = data_dir
script_dir = CDIR+'/XS_transformation_pHe'

sys.path.append( script_dir )



import lhcb
lhcb.write( data_dir,  res_dir )
