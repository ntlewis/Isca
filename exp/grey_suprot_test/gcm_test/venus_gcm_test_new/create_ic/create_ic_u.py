import analyse_functions_neil as af 
import netCDF4 as nc 
import numpy as np
import matplotlib.pyplot as plt





start =50
end = 50
exp = 'venus_col_log_check_newFlux_52H145'
ds = af.open_experiment(exp, start, end, 'atmos_monthly.nc', variables=['temp', 'ps'], cal=False)

# ds_lat = af.open_experiment('VENUS_T10L52_topsponge_full', 1, 1, 'atmos_monthly.nc', variables=['temp'], cal=False)
# lat = ds_lat.lat.values 
# np.save('T10_lat',lat)

lat = np.load('T21_lat.npy')
nlat = len(lat)
nlon = nlat*2

pfull = ds.pfull.values 
ps_cnst = 92e3
npres = len(pfull)

log_pres = -np.log(pfull/ps_cnst)
x = log_pres/log_pres[10]
xmax = x[0] 
umax = 100.
u1 = umax*x 
u2 = umax / (1-xmax) * (x-xmax)
# A = 100 * (B/(1-B)-1)
# C = 100*B / (1-B)
print(log_pres/log_pres[10], np.arange(len(log_pres)))
print(u1, u2, np.where(u1<u2, u1, u2))
uvert = np.where(u1<u2, u1, u2)
latmult = (np.abs(lat)-90)/(45-90)
latmult = np.where(latmult<1, latmult, 1)




print(latmult.shape, uvert.shape)

uprof = uvert[:,None] * latmult[None,:]

plt.contourf(lat, pfull, uprof)
plt.ylim([pfull[-1],pfull[0]])
plt.yscale('log')
plt.show()

#random numbers 
rng = np.random.default_rng(seed=42)
trn = (rng.random((npres,nlat,nlon)) - .5) / 1e3
urn = (rng.random((npres,nlat,nlon)) - .5) / 1e3
vrn = (rng.random((npres,nlat,nlon)) - .5) / 1e3
tsrn = (rng.random((nlat,nlon)) - .5) / 1e3
psrn = (rng.random((nlat,nlon)) - .5) / 1e3


print(ds)
ti  = ds.temp.values 
ps = ds.ps.values 
ti  = np.squeeze(ti[-1,:])
ps = np.squeeze(ps[-1])
# ti  = np.squeeze(ti[-1,:,:,0])
# ps = np.squeeze(ps[-1,:,0])

t  = np.ones((npres,nlat,nlon)) * ti[:,None,None] + trn#[:,:,None]
psurf = np.ones((nlat,nlon)) * ps[None,None] + psrn#[:,None]
tsurf = np.ones((nlat,nlon)) * np.squeeze(ti[-1])[None,None] + tsrn 

print(np.shape(tsurf), np.shape(np.squeeze(ti[-1])))

u  = np.ones_like(t)*uprof[:,:,None] + urn
v  = np.zeros_like(t) + vrn 
sh = np.zeros_like(t)

f = nc.Dataset(exp+'_icwithu_T21.nc','w',format='NETCDF3_CLASSIC')

f.createDimension('lon'  , nlon) 
f.createDimension('lat'  , nlat)
f.createDimension('pfull', npres)

longitude = f.createVariable('lon'   , 'f4', 'lon')
latitude  = f.createVariable('lat'   , 'f4', 'lat')
pfull     = f.createVariable('pfull' , 'f4', 'pfull')
temp      = f.createVariable('temp'  , 'f4', ('pfull', 'lat', 'lon'))
ucomp     = f.createVariable('ucomp' , 'f4', ('pfull', 'lat', 'lon'))
vcomp     = f.createVariable('vcomp' , 'f4', ('pfull', 'lat', 'lon'))
sphum     = f.createVariable('sphum' , 'f4', ('pfull', 'lat', 'lon'))
ps        = f.createVariable('ps'    , 'f4', ('lat', 'lon'))
t_surf    = f.createVariable('t_surf', 'f4', ('lat', 'lon'))

print(np.shape(t))

longitude[:] = np.arange(0,360,360/nlon)
latitude[:]  = lat
pfull[:]     = ds.pfull.values 
temp[:,:,:]  = t
ucomp[:,:,:] = u 
vcomp[:,:,:] = v 
sphum[:,:,:] = sh
ps[:,:]      = psurf
t_surf[:,:]  = tsurf

f.close()
