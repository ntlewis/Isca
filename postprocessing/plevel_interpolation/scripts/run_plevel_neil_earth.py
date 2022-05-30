from netCDF4 import Dataset  
from plevel_fn import plevel_call
import sys
import os
import time
import pdb
import subprocess


base_dir='/network/group/aopp/planetary/PLR016_LEWIS_VENUS/isca_data/'
out_dir=base_dir
#exp_name_list = ['omg_1_1_T85_dailyfinal2','omg_1_4_T85_dailyfinal2','omg_1_8_T85_dailyfinal2','omg_1_16_T85_dailyfinal2','omg_1_32_T85_dailyfinal2','omg_1_64_T85_dailyfinal2','omg_1_256_T85_dailyfinal2']
#exp_name_list = ['1','4','8','16','32','64','256']
#exp_extra = lambda x: 'omg_1_'+x+'_T85L30_axissym'
#exp_name_list = ['realistic_continents_fixed_sst_test_experiment']
# exp_name_list = ['omg_8_1_T127L80_FINAL','omg_4_1_T127L80_FINAL','omg_2_1_T127L80_FINAL','omg_1_1_T42L30_FINAL','omg_1_2_T42L30_FINAL','omg_1_4_T42L30_FINAL','omg_1_8_T42L30_FINAL','omg_1_16_T42L30_FINAL','omg_1_32_T42L30_FINAL','omg_1_64_T42L30_FINAL','omg_1_128_T42L30_FINAL']
# start_files=[33, 33, 33, 33, 33, 33, 67, 85, 85, 85, 85]
# end_files=  [36, 36, 36, 36, 36, 36, 72, 96, 96, 96, 96]

# exp_name_list2 = ['omg_4_1_T127L80_FINAL','omg_2_1_T127L80_FINAL','omg_1_1_T42L30_FINAL','omg_1_2_T42L30_FINAL','omg_1_4_T42L30_FINAL','omg_1_8_T42L30_FINAL','omg_1_16_T42L30_FINAL','omg_1_32_T42L30_FINAL','omg_1_64_T42L30_FINAL']
# start_files2=[ 1, 1, 1, 1, 1, 1, 1, 1, 1]
# end_files2=  [18,18,18,18,18,18,18,18,18]


# exp_name_list = ['omg_4_1_T127L80ax_FINAL','omg_2_1_T127L80ax_FINAL','omg_1_1_T127L80ax_FINAL','omg_1_2_T127L80ax_FINAL','omg_1_4_T127L80ax_FINAL','omg_1_8_T127L80ax_FINAL','omg_1_16_T127L80ax_FINAL','omg_1_32_T127L80ax_FINAL','omg_1_64_T127L80ax_FINAL','omg_1_128_T127L80ax_FINAL','omg_1_256_T127L80ax_FINAL']
# start_files=[ 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25]
# end_files=  [36,36,36,36,36,36,36,36,36, 36, 36]


exp_name_list = ['sbm_realistic']#,'omg_1_4_T127L80_FINAL','omg_1_64_T127L80_FINAL']#,'omg_1_16_T42L80ax_FINAL']
start_files = [10]#,  1,  1,  1]#,1]
end_files   = [10]#, 18, 18, 18]#,24]

###!!!!!! THIS HAS BEEN CHANGED TO -X
mask_below_surface_set=' ' #Default is to mask values that lie below the surface pressure when interpolated. For some applications, e.g. Tom Clemo's / Mark Baldwin's stratosphere index, you want to have values interpolated below ground, i.e. as if the ground wasn't there. To use this option, this value should be set to '-x '. 


# plevs =' -0 -p "200 800 1400 2200 3200 4500 5900 7400 9100 10800 12600 14500 16500 18600 20800 23200 25700 28600 31600 35000 38700 42800 47300 52200 57700 63800 70500 77900 86100 95200"'
# var_names='ucomp vcomp height temp omega udt_rdamp dt_ug_damp'
# file_suffix='_interp'

#plevs =' -0 -p "2 16 51 138 324 676 1266 2162 3407 5014 6957 9185 11627 14210 16864 19534 22181 24783 27331 29830 32290 34731 37173 39637 42147 44725 47391 50164 53061 56100 59295 62661 66211 66958 73915 78095 82510 87175 92104 97312"'
var_names='ucomp vcomp temp height sphum'
file_suffix='_interp'

plevs =' -p "2 16 51 138 324 676 1266 2162 3407 5014 6957 9185 11627 14210 16864 19534 22181 24783 27331 29830 32290 34731 37173 39637 42147 44725 47391 50164 53061 56100 59295 62661 66211 66958 73915 78095 82510 87175 92104 97312"'

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