#! /usr/bin/env python

import argparse
from PlotFunctions import plot2D


parser = argparse.ArgumentParser(description='Michael Korsmeier - Plotting script for cross section tables, in the form requiered by AntiDeuteronM')

parser.add_argument('--file',  help='Filename of table.', action='store', dest='file', type=str, default='')
parser.add_argument('--file2', help='Filename of second table (in case of ratio)', action='store', dest='file2', type=str, default='')
parser.add_argument('--o',     help='Filename of result file', action='store', dest='o', type=str, default='result.png')

parser.add_argument('--Tmin_proj',  help='Plotting range.', action='store', dest='Tmin_proj', type=float, default=1e0)
parser.add_argument('--Tmax_proj',  help='Plotting range.', action='store', dest='Tmax_proj', type=float, default=1e7)
parser.add_argument('--Tmin_prod',  help='Plotting range.', action='store', dest='Tmin_prod', type=float, default=1e-1)
parser.add_argument('--Tmax_prod',  help='Plotting range.', action='store', dest='Tmax_prod', type=float, default=1e7)

parser.add_argument('--Zmin',  help='Plotting range.', action='store', dest='Zmin', type=float, default=-1)
parser.add_argument('--Zmax',  help='Plotting range.', action='store', dest='Zmax', type=float, default=-1)

parser.add_argument('--cmap',  help='Color map.', action='store', dest='cmap', type=str, default='d')


args = parser.parse_args()

file   	    = args.file
file2 	    = args.file2
Tmin_proj   = args.Tmin_proj
Tmax_proj   = args.Tmax_proj
Tmin_prod   = args.Tmin_prod
Tmax_prod   = args.Tmax_prod
Zmin        = args.Zmin
Zmax        = args.Zmax
cmap_str    = args.cmap
resfile     = args.o


plot2D(file, file2, Tmin_proj, Tmax_proj, Tmin_prod, Tmax_prod, Zmin, Zmax, cmap_str, resfile)
