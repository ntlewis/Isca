import analyse_functions_neil as af 
import matplotlib.pyplot as plt 
import numpy as np

exp1 = 'column_test_rrtm_ozonefile_manylat_ob23'
exp2 = 'column_test_rrtm_nofile_manylat_ob23'
exp3 = 'column_test_rrtm_nofile_manylat'
exp4 = 'column_test_rrtm_ozonefile_manylat'

ds1 = af.open_experiment(exp1, 24, 32, 'atmos_monthly.nc')
ds2 = af.open_experiment(exp2, 24, 36, 'atmos_monthly.nc')
ds3 = af.open_experiment(exp3, 12, 22, 'atmos_monthly.nc')
ds4 = af.open_experiment(exp4, 12, 22, 'atmos_monthly.nc')

af.global_average_lat_lon(ds1, 'temp')
af.global_average_lat_lon(ds2, 'temp')
af.global_average_lat_lon(ds3, 'temp')
af.global_average_lat_lon(ds4, 'temp')



plt.figure(1)
for n, lat in enumerate(ds4.lat.values):
    plt.plot(np.squeeze(ds4.temp.mean('time').values[:,n]), ds4.pfull.values, color='b', alpha=0.3)
plt.plot(np.squeeze(ds4.temp_area_av.mean('time').values), ds4.pfull.values, color='b', linewidth=5)
for n, lat in enumerate(ds1.lat.values):
    plt.plot(np.squeeze(ds1.temp.mean('time').values[:,n]), ds1.pfull.values, color='g', alpha=0.3)
plt.plot(np.squeeze(ds1.temp_area_av.mean('time').values), ds1.pfull.values, color='g', linewidth=5)
for n, lat in enumerate(ds4.lat.values):
    plt.plot(np.squeeze(ds3.temp.mean('time').values[:,n]), ds3.pfull.values, color='r', alpha=0.3)
plt.plot(np.squeeze(ds3.temp_area_av.mean('time').values), ds3.pfull.values, color='r', linewidth=5)
for n, lat in enumerate(ds4.lat.values):
    plt.plot(np.squeeze(ds2.temp.mean('time').values[:,n]), ds2.pfull.values, color='k', alpha=0.3)
plt.plot(np.squeeze(ds2.temp_area_av.mean('time').values), ds2.pfull.values, color='k', linewidth=5)
plt.legend(loc='best')
plt.ylim([ds4.pfull.values.max(), ds4.pfull.values.min()])
plt.yscale('log')



plt.show()
