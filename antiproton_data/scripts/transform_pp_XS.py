#! /usr/bin/env python

import os
import sys
import glob
import math
import numpy as np



CDIR = os.path.dirname(os.path.realpath(__file__))

data_dir = CDIR+'/../data/pp'
res_dir  = data_dir
script_dir = CDIR+'/XS_transformation_pp'

sys.path.append( script_dir )

print( data_dir )

import na61
na61.write( data_dir+'/na61',               res_dir )

import allaby
allaby.write( data_dir+'/allaby',           res_dir )

import johnson
johnson.write( data_dir+'/johnson',         res_dir )

import antreasyan
antreasyan.write( data_dir+'/antreasyan',   res_dir )

import guettler
guettler.write( data_dir+'/guettler',       res_dir )

import capiluppi
capiluppi.write( data_dir+'/capiluppi',     res_dir )

import dekkers
dekkers.write( data_dir+'/dekkers',         res_dir )

import na49
na49.write( data_dir+'/na49',               res_dir )

import brahms
brahms.write( data_dir+'/brahms',           res_dir )

import phenix
phenix.write( data_dir+'/phenix',           res_dir )



