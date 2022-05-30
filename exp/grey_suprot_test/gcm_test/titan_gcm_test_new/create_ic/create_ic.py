import analyse_functions_neil as af 
import netCDF4 as nc 
import numpy as np






start =30
end = 30
exp = 'titan_col_log_check_newFlux_52H145_multicol'
ds = af.open_experiment(exp, start, end, 'atmos_monthly.nc', variables=['temp', 'ps'], cal=False)

# ds_lat = af.open_experiment('TITAN_T21L50', 1, 1, 'atmos_monthly.nc', variables=['temp'], cal=False)
# lat = ds_lat.lat.values 
# np.save('T21_lat',lat)

lat = np.load('T21_lat.npy')
nlat = len(lat)
nlon = nlat*2

pfull = ds.pfull.values 
npres = len(pfull)


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

t  = np.ones((npres,nlat,nlon)) * ti[:,:,None] + trn#[:,:,None]
psurf = np.ones((nlat,nlon)) * ps[:,None] + psrn#[:,None]
tsurf = np.ones((nlat,nlon)) * np.squeeze(ti[-1])[:,None] + tsrn 

print(np.shape(tsurf), np.shape(np.squeeze(ti[-1])))

u  = np.zeros_like(t) + urn
v  = np.zeros_like(t) + vrn 
sh = np.zeros_like(t)

f = nc.Dataset(exp+'_ic_T21lat.nc','w',format='NETCDF3_CLASSIC')

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
