#!/usr/bin/env python

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# wrf_mslp.py
# 
# Author(s): Jon Seibert
#
# Copyright (C) [2023]-[2033] Jon Seibert, Department of Meteorology & Atmospheric Science, The Pennsylvania State University, USA - All Rights Reserved 
#
# Last Modified: 17 January 2023
#
# Purpose: To calculate the Mean Sea Level Pressure (MSLP) of a WRF output file, and output the result as a new NetCDF file
#
# Usage: 
#			
# Usage examples:

# Inputs: * [One NetCDF file]
#
# Outputs: *.nc [One NetCDF file]
#			
# Package requirements: netCDF4, numpy, os, scipy, sys
#
# Notes: # THIS MIGHT ACTUALLY JUST BE GEOPOTENTIAL- IF SO CONVERT TO GPH BY DIVIDING GP BY g0 (9.81 etc)
# 		- 
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


## Imports
import netCDF4
from netCDF4 import Dataset
import numpy as np
import os
from os import listdir
import scipy.io.netcdf as nc
import sys
import wrf

## SETTINGS

input_filename = 'wrf_enkf_output_d02_020'
input_filepath = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/input/wrf_data/2022/1200_test'
output_filepath = 'C:/Users/Jon/Documents/Actual Documents/PSU/IMPACTS/Code/output/mslp'
output_filename = 'wrf_enkf_output_d02_020_mslp.nc'

lon_name = 'XLONG'
lat_name = 'XLAT'
data_name_refl = 'REFL_10CM' # Variable being analyzed
data_name_ppt = 'T' # PERTURBATION POTENTIAL TEMPERATURE
data_name_pres = 'PB' # Base Pressure
data_name_pres_2 = 'P' # Perturbation pressure
data_name_gph = 'PH' # Base geopotential height
data_name_gph_2 = 'PHB' # Perturbation geopotential height
data_name_q = 'QVAPOR' # Water vapor mixing ratio


## MAIN

# Create full target filepath
input_filepath_full = input_filepath + input_filename
output_filepath_full = output_filepath + output_filename

# Read in
raw_nc = Dataset(input_filepath_full)

# Extract variables
lon = raw_nc.variables[lon_name][:]
lat = raw_nc.variables[lat_name][:]
ppt = raw_nc.variables[data_name_ppt][:]
pressure_base = raw_nc.variables[data_name_pres][:]
pressure_pert = raw_nc.variables[data_name_pres_2][:]
gph_base = raw_nc.variables[data_name_gph][:]
gph_pert = raw_nc.variables[data_name_gph_2][:]
gph_pert = raw_nc.variables[data_name_gph_2][:]
mixing_ratio = raw_nc.variables[data_name_q][:]

# Compute input variables

data_pressure = pressure_base + pressure_pert # Pressure
data_gph = gph_base + gph_pert # Geopotential height  # THIS MIGHT ACTUALLY JUST BE GEOPOTENTIAL- IF SO CONVERT TO GPH BY DIVIDING GP BY g0 (9.81 etc)

pt = ppt + 300
data_t = np.multiply(pt,np.power((data_pressure/100000),0.287) # Temperature

data_q = mixing_ratio; # Water vapor mixing ratio

# Compute MSLP

mslp = wrf.slp(data_gph,data_t,data_pressure,data_q,False,'mb')

# Write out MSLP as NetCDF file with one variable

ncfile = Dataset(output_filepath_full,mode='w',format='NETCDF4_CLASSIC')
lat_dim = ncfile.createDimension('XLAT',399)
lon_dim = ncfile.createDimension('XLONG',399)
time_dim = ncfile.createDimension('TIME',1)

# Create variables

lat = ncfile.createVariable('XLAT', np.float32, ('lat',))
lat.units = 'degrees_north'
lat.long_name = 'latitude'
lon = ncfile.createVariable('XLONG', np.float32, ('lon',))
lon.units = 'degrees_west'
lon.long_name = 'longitude'
mslp_var = ncfile.createVariable('MSLP', np.float32, ('lon','lat'))
mslp_var.units = 'Pa'
mslp_var.long_name = 'Mean Sea Level Pressure (MSLP)'

#lat = ncfile.createVariable('lat', np.float32, ('lat',))
#lat.units = 'degrees_north'
#lat.long_name = 'latitude'

# Apparently I can compute this just with the geopotential height, because gph is just the height at which the pressure is flat or some shit idk google it