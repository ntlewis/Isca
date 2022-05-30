import analyse_functions_neil as af 
import matplotlib.pyplot as plt 
import numpy as np

# exp1 = 'column_test_wts1'
# exp2 = 'column_test_wts100'
#exp3 = 'column_test_rrtm_noozone'
exp4 = 'column_test_exp_remspec'

# ds1 = af.open_experiment(exp1, 37, 48, 'atmos_monthly.nc')
# ds2 = af.open_experiment(exp2, 37, 48, 'atmos_monthly.nc')
#ds3 = af.open_experiment(exp3, 10, 10, 'atmos_monthly.nc', variables=['temp','coszen'], cal=False)
ds4 = af.open_experiment(exp4, 12, 12, 'atmos_monthly.nc', variables=['temp'], cal=False)

# plt.figure(1)
# plt.plot(np.squeeze(ds1.temp.mean('time').values), ds1.pfull.values, label='wts1')
# plt.plot(np.squeeze(ds2.temp.mean('time').values), ds2.pfull.values, label='wts100')
# plt.legend(loc='best')
# plt.ylim([ds1.pfull.values.max(), ds1.pfull.values.min()])
# plt.yscale('log')

temp = np.array([291.4, 288.1, 284.9, 281.7, 278.4, 275.2, 271.9, 268.7, 265.4, 262.2, 258.9, 255.7, 252.4, 249.2, 245.9, 242.7, 239.5, 236.2, 233.0, 229.7, 226.5, 223.3, 220.0, 216.8, 216.6, 216.6, 216.6, 216.6])
patmos = np.array([107477, 101325, 95461, 89876, 84559, 79501, 74691, 70121, 65780, 61660, 57752, 54048, 50539, 47217, 44075, 41105, 38299, 35651, 33154, 30800, 28584, 26499, 24540, 22699, 20984, 19399, 17933, 16579])


plt.figure(2, figsize=(7,5))
plt.title('Grey Experiment', fontsize=15)
plt.plot(np.squeeze(ds4.temp.mean('time').values)[1:], ds4.pfull.values[1:], linewidth=3)
plt.plot(temp, patmos/100., label='US STD ATMOS', linewidth=3, color='grey', alpha=0.5)
plt.ylim([ds4.pfull.values.max(), ds4.pfull.values[1:].min()])
plt.yscale('log')
plt.ylabel('Pressure (hPa)', fontsize=13)
plt.xlabel('Temperature (K)', fontsize=13)



plt.show()
