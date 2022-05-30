from netCDF4 import Dataset  
from plevel_fn import plevel_call
import sys
import os
import time
import pdb
import subprocess

start_time=time.time()
base_dir='/network/group/aopp/planetary/PLR018_LEWIS_VENUS2/isca_data/'
out_dir=base_dir
#exp_name_list = ['omg_1_1_T85_dailyfinal2','omg_1_4_T85_dailyfinal2','omg_1_8_T85_dailyfinal2','omg_1_16_T85_dailyfinal2','omg_1_32_T85_dailyfinal2','omg_1_64_T85_dailyfinal2','omg_1_256_T85_dailyfinal2']
#exp_name_list = ['1','4','8','16','32','64','256']
#exp_extra = lambda x: 'omg_1_'+x+'_T85L30_axissym'
#exp_name_list = ['realistic_continents_fixed_sst_test_experiment']
exp_name_list = ['omg_1_16_T170L80']
start_file=133
end_file=136
nfiles=(end_file-start_file)+1


mask_below_surface_set=' ' #Default is to mask values that lie below the surface pressure when interpolated. For some applications, e.g. Tom Clemo's / Mark Baldwin's stratosphere index, you want to have values interpolated below ground, i.e. as if the ground wasn't there. To use this option, this value should be set to '-x '. 


# plevs =' -0 -p "200 800 1400 2200 3200 4500 5900 7400 9100 10800 12600 14500 16500 18600 20800 23200 25700 28600 31600 35000 38700 42800 47300 52200 57700 63800 70500 77900 86100 95200"'
# var_names='ucomp vcomp height temp omega udt_rdamp dt_ug_damp'
# file_suffix='_interp'

#plevs =' -0 -p "2 16 51 138 324 676 1266 2162 3407 5014 6957 9185 11627 14210 16864 19534 22181 24783 27331 29830 32290 34731 37173 39637 42147 44725 47391 50164 53061 56100 59295 62661 66211 66958 73915 78095 82510 87175 92104 97312"'
var_names='ucomp vcomp height temp omega'
file_suffix='_interp'

plevs =' -0 -p "124 391 515 666 846 1056 1298 1574 1883 2224 2599 3004 3439 3902 4390 4903 5437 5991 6563 7151 7753 8369 8996 9635 10284 10943 11613 12292 12982 13683 14396 15121 15860 16612 17381 18166 18970 19793 20638 21505 22396 23314 24259 25233 26238 27277 28350 29459 30608 31796 33027 34302 35624 36994 38415 39889 41418 43004 44651 46359 48132 49973 51883 53866 55925 58063 60281 62585 64977 67460 70037 72714 75492 78377 81372 84481 87710 91061 94541 98153"'


for exp_name in exp_name_list:
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
        nc_file_in = base_dir+'/'+exp_name+'/run'+number_prefix+str(n+start_file)+'/atmos_out.nc'
        nc_file_out = out_dir+'/'+exp_name+'/run'+number_prefix+str(n+start_file)+'/atmos_out'+file_suffix+'.nc'

        if not os.path.isfile(nc_file_out):
            plevel_call(nc_file_in,nc_file_out, var_names = var_names, p_levels = plevs, mask_below_surface_option=mask_below_surface_set)             



print('execution time', time.time()-start_time)



