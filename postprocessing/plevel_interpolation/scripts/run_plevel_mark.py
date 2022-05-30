from netCDF4 import Dataset  
from plevel_fn import plevel_call
import sys
import os
import time
import pdb
import subprocess

start_time=time.time()
#base_dir='/network/group/aopp/planetary/PLR016_LEWIS_VENUS/isca_data/omg_1_1_T85_dailyfinal2/run0001/'
base_dir='/network/group/aopp/planetary/RTP002_HAMMOND_TL-DYN/exofms-cs-output/neil-alpha-10day-cold/neil-alpha-10day-cold_2000/'
#out_dir='/network/group/aopp/planetary/PLR018_LEWIS_VENUS2/'

exp_name_list = ['']

mask_below_surface_set=' ' #Default is to mask values that lie below the surface pressure when interpolated. For some applications, e.g. Tom Clemo's / Mark Baldwin's stratosphere index, you want to have values interpolated below ground, i.e. as if the ground wasn't there. To use this option, this value should be set to '-x '. 


#plevs =' -0 -p "200 800 1400 2200 3200 4500 5900 7400 9100 10800 12600 14500 16500 18600 20800 23200 25700 28600 31600 35000 38700 42800 47300 52200 57700 63800 70500 77900 86100 95200"'
#plevs =' -0 -p "8300 9500 10900 12500 14200 16200 18500 21100 24000 27200 30800 34700 39100 43900 49000 54600 60500 66400 72000 76900 81200"'
plevs =' -0 -p "70 150 270 440 670 940 1230 1550 1890 2270 2680 3120 3620 4180 4790 5470 6210 7030 7920 8900 9990 11210 12570 14110 15830 17760 19930 22360 25090 28150 31590 35420 39620 44180 49110 54390 59990 65810 71650 77140 81810 85560 88630 91260 93620 95740 97630 99250"'
var_names='ucomp vcomp height temp omega ppt'
file_suffix='_interp'



for exp_name in exp_name_list:

    #nc_file_in = base_dir+'/atmos_out.nc'
    #nc_file_out = out_dir+'/neil_atmos_daily_interp.nc'
    nc_file_in  = base_dir+'/atmos_daily.nc'
    nc_file_out = base_dir+'/atmos_daily_interp.nc'

    if not os.path.isfile(nc_file_out):
        plevel_call(nc_file_in,nc_file_out, var_names = var_names, p_levels = plevs, mask_below_surface_option=mask_below_surface_set)             



print('execution time', time.time()-start_time)



