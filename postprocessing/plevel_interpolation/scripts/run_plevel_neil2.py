from netCDF4 import Dataset  
from plevel_fn import plevel_call
import sys
import os
import time
import pdb
import subprocess


base_dir='/network/group/aopp/planetary/PLR018_LEWIS_VENUS2/isca_data/'
out_dir=base_dir

exp_name_list = ['base_exp_T42L20_om_1_1_a_1_1_ps_1_1_dsol_1_sc_1000_tau43_spg_dt360_dailyfin']#,'omg_1_4_T127L80_FINAL','omg_1_64_T127L80_FINAL']#,'omg_1_16_T42L80ax_FINAL']
start_files = [1]#,  1,  1,  1]#,1]
end_files   = [1]#, 18, 18, 18]#,24]

###!!!!!! THIS HAS BEEN CHANGED TO -X
mask_below_surface_set=' ' #Default is to mask values that lie below the surface pressure when interpolated. For some applications, e.g. Tom Clemo's / Mark Baldwin's stratosphere index, you want to have values interpolated below ground, i.e. as if the ground wasn't there. To use this option, this value should be set to '-x '. 


# plevs =' -0 -p "200 800 1400 2200 3200 4500 5900 7400 9100 10800 12600 14500 16500 18600 20800 23200 25700 28600 31600 35000 38700 42800 47300 52200 57700 63800 70500 77900 86100 95200"'
# var_names='ucomp vcomp height temp omega udt_rdamp dt_ug_damp'
# file_suffix='_interp'

#plevs =' -0 -p "2 16 51 138 324 676 1266 2162 3407 5014 6957 9185 11627 14210 16864 19534 22181 24783 27331 29830 32290 34731 37173 39637 42147 44725 47391 50164 53061 56100 59295 62661 66211 66958 73915 78095 82510 87175 92104 97312"'
var_names='ucomp vcomp omega temp height'
file_suffix='_interp'

#plevs =' -0 -p "67 147 271 446 671 938 1234 1552 1897 2272 2679 3127 3624 4179 4793 5471 6214 7030 7923 8904 9993 11213 12581 14117 15839 17772 19941 22374 25104 28167 31604 35436 39639 44203 49127 54411 60012 65828 71670 77150 81816 85569 88637 91261 93617 95746 97630 99254"'

plevs = ' -0 -p "1839 7358 12416 17440 22454 27462 32468 37472 42475 47478 52480 57482 62483 67485 72486 77587 82487 87488 92489 97489"'

start_time=time.time()
for j, exp_name in enumerate(exp_name_list):
    start_file = start_files[j]
    end_file = end_files[j]
    nfiles=(end_file-start_file)+1
    for n in range(nfiles):
        print(n+start_file)

        number_prefix=''

        if n+start_file < 1000:
            number_prefix='0'
        if n+start_file < 100:
            number_prefix='00'
        if n+start_file < 10:
            number_prefix='000'

        #nc_file_in = base_dir+'/'+exp_extra(exp_name)+'/run'+number_prefix+str(n+start_file)+'/atmos_out.nc'
        #nc_file_out = out_dir+'/'+exp_extra(exp_name)+'/run'+number_prefix+str(n+start_file)+'/atmos_out'+file_suffix+'.nc'
        nc_file_in = base_dir+'/'+exp_name+'/run'+number_prefix+str(n+start_file)+'/atmos_monthly.nc'
        nc_file_out = out_dir+'/'+exp_name+'/run'+number_prefix+str(n+start_file)+'/atmos_monthly'+file_suffix+'.nc'

        if not os.path.isfile(nc_file_out):
            plevel_call(nc_file_in,nc_file_out, var_names = var_names, p_levels = plevs, mask_below_surface_option=mask_below_surface_set)             



print('execution time', time.time()-start_time)

# start_time=time.time()
# for j, exp_name in enumerate(exp_name_list2):
#     start_file = start_files2[j]
#     end_file = end_files2[j]
#     nfiles=(end_file-start_file)+1
#     for n in range(nfiles):
#         print(n+start_file)

#         number_prefix=''

#         if n+start_file < 1000:
#             number_prefix='0'
#         if n+start_file < 100:
#             number_prefix='00'
#         if n+start_file < 10:
#             number_prefix='000'

#         #nc_file_in = base_dir+'/'+exp_extra(exp_name)+'/run'+number_prefix+str(n+start_file)+'/atmos_out.nc'
#         #nc_file_out = out_dir+'/'+exp_extra(exp_name)+'/run'+number_prefix+str(n+start_file)+'/atmos_out'+file_suffix+'.nc'
#         nc_file_in = base_dir+'/'+exp_name+'/run'+number_prefix+str(n+start_file)+'/atmos_out.nc'
#         nc_file_out = out_dir+'/'+exp_name+'/run'+number_prefix+str(n+start_file)+'/atmos_out'+file_suffix+'.nc'

#         if not os.path.isfile(nc_file_out):
#             plevel_call(nc_file_in,nc_file_out, var_names = var_names, p_levels = plevs, mask_below_surface_option=mask_below_surface_set)       
# print('execution time2', time.time()-start_time)


# start_time=time.time()
# for j, exp_name in enumerate(exp_name_list3):
#     start_file = start_files3[j]
#     end_file = end_files3[j]
#     nfiles=(end_file-start_file)+1
#     for n in range(nfiles):
#         print(n+start_file)

#         number_prefix=''

#         if n+start_file < 1000:
#             number_prefix='0'
#         if n+start_file < 100:
#             number_prefix='00'
#         if n+start_file < 10:
#             number_prefix='000'

#         #nc_file_in = base_dir+'/'+exp_extra(exp_name)+'/run'+number_prefix+str(n+start_file)+'/atmos_out.nc'
#         #nc_file_out = out_dir+'/'+exp_extra(exp_name)+'/run'+number_prefix+str(n+start_file)+'/atmos_out'+file_suffix+'.nc'
#         nc_file_in = base_dir+'/'+exp_name+'/run'+number_prefix+str(n+start_file)+'/atmos_out.nc'
#         nc_file_out = out_dir+'/'+exp_name+'/run'+number_prefix+str(n+start_file)+'/atmos_out'+file_suffix+'.nc'

#         if not os.path.isfile(nc_file_out):
#             plevel_call(nc_file_in,nc_file_out, var_names = var_names, p_levels = plevs, mask_below_surface_option=mask_below_surface_set)       
# print('execution time3', time.time()-start_time)
