from netCDF4 import Dataset  
from plevel_fn import plevel_call
import sys
import os
import time
import pdb
import subprocess


base_dir='/network/group/aopp/planetary/PLR018_LEWIS_VENUS2/titan_data/'
out_dir=base_dir

exp_name_list = ['TAM_data_newz_nointerp']


###!!!!!! THIS HAS BEEN CHANGED TO -X
mask_below_surface_set=' ' #Default is to mask values that lie below the surface pressure when interpolated. For some applications, e.g. Tom Clemo's / Mark Baldwin's stratosphere index, you want to have values interpolated below ground, i.e. as if the ground wasn't there. To use this option, this value should be set to '-x '. 


var_names='ucomp vcomp omega temp zcomp'
file_suffix='_interp'

plevs =' -0 -p "1 2 4 11 25 55 115 232 447 826 1464 2491 4072 6405 9705 14183 20014 27310 36082 46226 57510 69591 82044 94406 106235 117159 126931 135477 142967"'

start_time=time.time()
for j, exp_name in enumerate(exp_name_list):

        nc_file_in = base_dir+'/'+exp_name+'.nc'
        nc_file_out = out_dir+'/'+exp_name+file_suffix+'.nc'

        if not os.path.isfile(nc_file_out):
            plevel_call(nc_file_in,nc_file_out, var_names = var_names, p_levels = plevs, mask_below_surface_option=mask_below_surface_set)             



print('execution time', time.time()-start_time)

