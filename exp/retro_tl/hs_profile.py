from __future__ import print_function
import numpy as np
from netCDF4 import Dataset, date2num
import pdb
import os
import xarray as xar

def output_to_file(data,lats,lons,latbs,lonbs,p_full,p_half,file_name,variable_name,number_dict):

	output_file = Dataset(file_name, 'w', format='NETCDF3_CLASSIC')

	is_thd=True


	lat = output_file.createDimension('lat', number_dict['nlat'])
	lon = output_file.createDimension('lon', number_dict['nlon'])

	latb = output_file.createDimension('latb', number_dict['nlatb'])
	lonb = output_file.createDimension('lonb', number_dict['nlonb'])

	if is_thd:
		pfull = output_file.createDimension('pfull', number_dict['npfull'])
		phalf = output_file.createDimension('phalf', number_dict['nphalf'])


	latitudes = output_file.createVariable('lat','d',('lat',))
	longitudes = output_file.createVariable('lon','d',('lon',))

	latitudebs = output_file.createVariable('latb','d',('latb',))
	longitudebs = output_file.createVariable('lonb','d',('lonb',))
	if is_thd:
		pfulls = output_file.createVariable('pfull','d',('pfull',))
		phalfs = output_file.createVariable('phalf','d',('phalf',))


	latitudes.units = 'degrees_N'.encode('utf-8')
	latitudes.cartesian_axis = 'Y'
	latitudes.edges = 'latb'
	latitudes.long_name = 'latitude'

	longitudes.units = 'degrees_E'.encode('utf-8')
	longitudes.cartesian_axis = 'X'
	longitudes.edges = 'lonb'
	longitudes.long_name = 'longitude'

	latitudebs.units = 'degrees_N'.encode('utf-8')
	latitudebs.cartesian_axis = 'Y'
	latitudebs.long_name = 'latitude edges'

	longitudebs.units = 'degrees_E'.encode('utf-8')
	longitudebs.cartesian_axis = 'X'
	longitudebs.long_name = 'longitude edges'

	if is_thd:
		pfulls.units = 'hPa'
		pfulls.cartesian_axis = 'Z'
		pfulls.positive = 'down'
		pfulls.long_name = 'full pressure level'

		phalfs.units = 'hPa'
		phalfs.cartesian_axis = 'Z'
		phalfs.positive = 'down'
		phalfs.long_name = 'half pressure level'

	if is_thd:
		output_array_netcdf = output_file.createVariable(variable_name,'f4',('pfull', 'lat','lon',))
	else:
		output_array_netcdf = output_file.createVariable(variable_name,'f4',('lat','lon',))

	latitudes[:] = lats
	longitudes[:] = lons

	latitudebs[:] = latbs
	longitudebs[:] = lonbs

	if is_thd:
		pfulls[:]     = p_full
		phalfs[:]     = p_half


	output_array_netcdf[:] = data

	output_file.close()

def hs_profile():

    # vertical coordinate
    etas = np.array([2.265e-6, 2.925e-6, 4.500e-6, 7.950e-6, 1.525e-5, 3.070e-5, 6.355e-5, 1.324e-4, 2.710e-4, 5.400e-4, 1.043e-3, 1.934e-3, 3.454e-3, 5.929e-3, 9.814e-3, 1.569e-2, 2.430e-2, 3.657e-2, 5.361e-2, 7.672e-2, 1.074e-1, 1.472e-1,1.979e-1, 2.610e-1, 3.373e-1, 4.273e-1, 5.299e-1, 6.420e-1, 7.577e-1, 8.679e-1, 9.602e-1])
    half_etas = np.concatenate(([0.],(etas[:-1] + etas[1:])/2.,[1.])) 
    print(half_etas, len(etas), len(half_etas))



    eq_pole_temp_grad = np.array([2.000, 2.370, 3.126, 4.439, 6.504, 9.599, 14.001, 19.807, 26.758, 34.295, 41.536, 47.342, 50.757, 51.175, 48.534, 43.315, 36.373, 28.746, 21.370, 14.948, 9.840, 7.524, 7.167, 6.762, 6.719, 6.678, 6.634, 6.410, 6.176, 6.136, 6.101])

    poly_coeffs = np.array([733.4861, 301.2670, 47.2002, -5.0230, -2.3257, -0.1896])

    lats = np.linspace(-89,89,101)


    lons = [0.]
    lonb_temp = [0., 360.]

    latb_temp=np.zeros(lats.shape[0]+1)

    for tick in np.arange(1,lats.shape[0]):
        latb_temp[tick]=(lats[tick-1]+lats[tick])/2.

    latb_temp[0]=-90.
    latb_temp[-1]=90.

    latbs=latb_temp
    lonbs=lonb_temp

    lats[0]=-90.+0.01
    lats[-1]=90.-0.01


    nlon=1#lons.shape[0]
    nlat=lats.shape[0]

    nlonb=len(lonbs)
    nlatb=latbs.shape[0]

    def poly_sum(eta):
    
        T_sum = 0
        for i in range(0,6):
            T_sum = T_sum + poly_coeffs[i] * (np.log10(eta))**i

        return T_sum

    T_eq = np.zeros([len(etas), len(lats), 1])
    T_res = np.zeros([len(etas), len(lats), 1])
    for j, lat in enumerate(lats):
        for k, eta in enumerate(etas):
            T_res[k, j, 0] = eq_pole_temp_grad[k] * (np.cos(np.deg2rad(lat)) - np.pi/4.)
            T_eq[k, j, 0] = poly_sum(eta) + T_res[k, j, 0]


    p_full = etas * 9.2e6 / 100.
    p_half = half_etas * 9.2e6 / 100.
    #print(p_full)

    time_arr = [0.]
    


    #DO FLIPPING

    #p_full=p_full[::-1]
    #p_half=p_half[::-1]

    #sh_new=sh_new[:,::-1,:]


    #Find grid and time numbers

    nlon=len(lons)
    nlat=len(lats)

    nlonb=len(lonbs)
    nlatb=len(latbs)

    npfull=len(p_full)
    nphalf=len(p_half)

    ntime=len(time_arr)


    #Output it to a netcdf file. 
    file_name='venus_relaxation_profile.nc'
    variable_name='venus_relaxation_profile'

    number_dict={}
    number_dict['nlat']=nlat
    number_dict['nlon']=nlon
    number_dict['nlatb']=nlatb
    number_dict['nlonb']=nlonb
    number_dict['npfull']=npfull
    number_dict['nphalf']=nphalf

    output_to_file(T_eq,lats,lons,latbs,lonbs,p_full,p_half, file_name,variable_name,number_dict)
    print('here')
if __name__ == "__main__":
    
    hs_profile()
