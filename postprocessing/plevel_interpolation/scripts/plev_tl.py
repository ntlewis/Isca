from netCDF4 import Dataset  
from plevel_fn import plevel_call
import sys
import os
import time
import pdb
import subprocess


base_dir='/network/group/aopp/planetary/PLR018_LEWIS_VENUS2/isca_data/tl_1_1rad_1_4rt2om_1_1tau_rad_diff10/run0030'
out_dir=base_dir

exp_name_list = ['atmos_monthly']


###!!!!!! THIS HAS BEEN CHANGED TO -X
mask_below_surface_set=' ' #Default is to mask values that lie below the surface pressure when interpolated. For some applications, e.g. Tom Clemo's / Mark Baldwin's stratosphere index, you want to have values interpolated below ground, i.e. as if the ground wasn't there. To use this option, this value should be set to '-x '. 


var_names='ucomp vcomp omega temp'
file_suffix='_interp'

plevs =' -0 -p "1839 7357 12415 17440 22453 27462 32467 37472 42475 47478 52480 57481 62483 67484 72485 77486 82487 87488 92488 97489"'

start_time=time.time()
for j, exp_name in enumerate(exp_name_list):

        nc_file_in = base_dir+'/'+exp_name+'.nc'
        nc_file_out = out_dir+'/'+exp_name+file_suffix+'.nc'

        if not os.path.isfile(nc_file_out):
            plevel_call(nc_file_in,nc_file_out, var_names = var_names, p_levels = plevs, mask_below_surface_option=mask_below_surface_set)             



print('execution time', time.time()-start_time)

