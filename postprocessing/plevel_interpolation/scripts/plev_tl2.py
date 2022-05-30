from netCDF4 import Dataset  
from plevel_fn import plevel_call
import sys
import os
import time
import pdb
import subprocess


base_dir='/network/group/aopp/planetary/PLR018_LEWIS_VENUS2/isca_data/F_greytl_1_1rad_'
out_dir=base_dir

f_name = 'atmos_monthly'
exp_name_list = [#'1_1om', '1_1po14w2om', '1_1rt2om',  '1_1po34w2om',
         #'1_2om','1_2po14w2om', '1_2rt2om',  '1_2po34w2om', 
         #'1_4om', '1_4po14w2om', '1_4rt2om',  '1_4po34w2om', 
         '1_16om',# '1_8po14w2om', '1_8rt2om',  '1_8po34w2om', 
         #'1_16om','1_16po14w2om','1_16rt2om', '1_16po34w2om', 
         #'1_32om','1_32po14w2om','1_32rt2om', '1_32po34w2om', 
         #'1_64om'
         ]
#run_name_list = ['run0020','run0021','run0022','run0023','run0024','run0025',
#                 'run0026','run0027',
#                 'run0028',
#                 'run0029', 'run0030'
#                ]
run_name_list = ['run0031','run0032','run0033','run0034','run0035','run0036']

###!!!!!! THIS HAS BEEN CHANGED TO -X
mask_below_surface_set=' ' #Default is to mask values that lie below the surface pressure when interpolated. For some applications, e.g. Tom Clemo's / Mark Baldwin's stratosphere index, you want to have values interpolated below ground, i.e. as if the ground wasn't there. To use this option, this value should be set to '-x '. 


var_names='ucomp vcomp omega temp height'
file_suffix='_interp'

plevs =' -0 -p "913 1554 2536 3972 5985 8690 12179 16508 21679 27631 34242 41335 48691 56071 63236 69970 76097 81493 86090 89878 92898 95237 97015 98377 99487"'

start_time=time.time()
for j, exp_name in enumerate(exp_name_list):
        for i, run_name in enumerate(run_name_list):

                nc_file_in = base_dir+exp_name+'/'+run_name+'/'+f_name+'.nc'
                nc_file_out = out_dir+exp_name+'/'+run_name+'/'+f_name+file_suffix+'.nc'

                if not os.path.isfile(nc_file_out):
                        plevel_call(nc_file_in,nc_file_out,
                                    var_names = var_names, p_levels = plevs,
                                    mask_below_surface_option=mask_below_surface_set)             



print('execution time', time.time()-start_time)

