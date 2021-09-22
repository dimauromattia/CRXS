#! /usr/bin/env python

import os
import sys
import glob
import math
import numpy as np


CRACS = os.environ['CRACS']
data_dir    = CRACS+'/data/CS_measurements/pBe'
res_dir     = data_dir
script_dir  = CRACS+'/scripts/CS_transformation_pBe'

sys.path.append( script_dir )


import veronin
veronin.write( data_dir+'/veronin',  res_dir )
