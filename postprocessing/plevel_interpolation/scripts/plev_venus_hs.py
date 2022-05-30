from netCDF4 import Dataset  
from plevel_fn import plevel_call
import sys
import os
import time
import pdb
import subprocess


base_dir='/network/group/aopp/planetary/PLR018_LEWIS_VENUS2/isca_data/venus_test_daily_inst/run0001'
out_dir=base_dir

exp_name_list = ['atmos_out']


###!!!!!! THIS HAS BEEN CHANGED TO -X
mask_below_surface_set=' ' #Default is to mask values that lie below the surface pressure when interpolated. For some applications, e.g. Tom Clemo's / Mark Baldwin's stratosphere index, you want to have values interpolated below ground, i.e. as if the ground wasn't there. To use this option, this value should be set to '-x '. 


var_names='ucomp vcomp omega temp'
file_suffix='_interp'

plevs =' -0 -p "8 28 45 80 156 315 653 1350 2739 5408 10321 18968 33552 57163 93971 149394 230305 345216 504350 719685 1004691 1374048 1842998 2424487 3126866 3951858 4888612 5906931 6951716 7938955 8801669"'

start_time=time.time()
for j, exp_name in enumerate(exp_name_list):

        nc_file_in = base_dir+'/'+exp_name+'.nc'
        nc_file_out = out_dir+'/'+exp_name+file_suffix+'.nc'

        if not os.path.isfile(nc_file_out):
            plevel_call(nc_file_in,nc_file_out, var_names = var_names, p_levels = plevs, mask_below_surface_option=mask_below_surface_set)             



print('execution time', time.time()-start_time)

