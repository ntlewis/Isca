from netCDF4 import Dataset  
from plevel_fn import plevel_call
import sys
import os
import time
import pdb
import subprocess


base_dir='/network/group/aopp/planetary/PLR018_LEWIS_VENUS2/titan_data/new_titan/'
out_dir=base_dir

exp_name_list = ['atmos_nointerp']


###!!!!!! THIS HAS BEEN CHANGED TO -X
mask_below_surface_set=' ' #Default is to mask values that lie below the surface pressure when interpolated. For some applications, e.g. Tom Clemo's / Mark Baldwin's stratosphere index, you want to have values interpolated below ground, i.e. as if the ground wasn't there. To use this option, this value should be set to '-x '. 


var_names='ucomp vcomp omega temp zg'
file_suffix='_interp'

plevs =' -0 -p "15352  18824  22810  27319  32350  37886  43896  50330  57127  64209  71492  78880  86276  93582 100707 107564 114083 120207 125898 131144 135954 140374 144495"'

start_time=time.time()
for j, exp_name in enumerate(exp_name_list):

        nc_file_in = base_dir+'/'+exp_name+'.nc'
        nc_file_out = out_dir+'/'+exp_name+file_suffix+'.nc'

        if not os.path.isfile(nc_file_out):
            plevel_call(nc_file_in,nc_file_out, var_names = var_names, p_levels = plevs, mask_below_surface_option=mask_below_surface_set)             



print('execution time', time.time()-start_time)

