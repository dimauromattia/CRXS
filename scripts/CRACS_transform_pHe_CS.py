#! /usr/bin/env python

import os
import sys
import glob
import math
import numpy as np


CRACS = os.environ['CRACS']
data_dir = CRACS+'/data/CS_measurements/pHe'
res_dir  = data_dir
script_dir = CRACS+'/scripts/CS_transformation_pHe'

sys.path.append( script_dir )


import lhcb
lhcb.write( data_dir,  res_dir )
