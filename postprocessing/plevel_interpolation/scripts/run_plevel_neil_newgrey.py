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
#exp_name_list = ['TITAN_T21L52_olddiff_highsponge_oldbl_diff3day_gust05_12hrly']
exp_name_list = ['TITAN_T21LLEB_dailyrestest60_dt120']
start_file=1
end_file=4
nfiles=(end_file-start_file)+1


mask_below_surface_set=' ' #Default is to mask values that lie below the surface pressure when interpolated. For some applications, e.g. Tom Clemo's / Mark Baldwin's stratosphere index, you want to have values interpolated below ground, i.e. as if the ground wasn't there. To use this option, this value should be set to '-x '. 


# plevs =' -0 -p "200 800 1400 2200 3200 4500 5900 7400 9100 10800 12600 14500 16500 18600 20800 23200 25700 28600 31600 35000 38700 42800 47300 52200 57700 63800 70500 77900 86100 95200"'
# var_names='ucomp vcomp height temp omega udt_rdamp dt_ug_damp'
# file_suffix='_interp'

#plevs =' -0 -p "2 16 51 138 324 676 1266 2162 3407 5014 6957 9185 11627 14210 16864 19534 22181 24783 27331 29830 32290 34731 37173 39637 42147 44725 47391 50164 53061 56100 59295 62661 66211 66958 73915 78095 82510 87175 92104 97312"'
var_names='ucomp vcomp omega ucomp_sq ucomp_vcomp ucomp_omega dt_ug_damp dt_ug_diffusion udt_rdamp'
file_suffix='_interp'

#plevs =' -0 -p "1 2 3 6 10 18 29 50 80 127 197 300 448 656 944 1332 1849 2521 3383 4467 5807 7438 9390 11691 14362 17419 20868 24708 28925 33500 38404 43601 49049 54701 60512 66433 72417 78425 84420 90374 96299 102087 107834 113516 119154 124776 130424 136149 142014"'

plevs =' -0 -p "1 2 4 5 7 11 15 21 30 41 58 81 114 159 221 309 432 603 842 1175 1641 2290 3192 4458 6226 8691 12131 16926 23631 32548 43338 55429 68434 81757 94722 106747 117269 126088 133072 138213 141642 143630 144590 144932"'

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
        nc_file_in = base_dir+'/'+exp_name+'/run'+number_prefix+str(n+start_file)+'/atmos_monthly.nc'
        nc_file_out = out_dir+'/'+exp_name+'/run'+number_prefix+str(n+start_file)+'/atmos_monthly'+file_suffix+'.nc'

        if not os.path.isfile(nc_file_out):
            plevel_call(nc_file_in,nc_file_out, var_names = var_names, p_levels = plevs, mask_below_surface_option=mask_below_surface_set)             



print('execution time', time.time()-start_time)



