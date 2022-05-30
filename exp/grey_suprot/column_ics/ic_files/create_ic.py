import analyse_functions_neil as af 
import netCDF4 as nc 
import numpy as np



# ds_lat = af.open_experiment('fr_T42_L50_dt30_damp2_5h_LAPsponge_0_4h_om_1_256', start, end, 'atmos_monthly.nc', variables=['temp'], cal=False)
# lat = ds_lat.lat.values 
# np.save('T42_lat',lat)

fracs = ['0','1_10','5_20','1_2','1']
pss = ['1','20']

for i in range(len(fracs)):
    for j in range(len(pss)):


        start =10
        end = 10
        exp = 'col_'+pss[j]+'bar_frac'+fracs[i]
        ds = af.open_experiment(exp, start, end, 'atmos_monthly.nc', variables=['temp', 'ps'], cal=False)

        lat = np.load('T42_lat.npy')

        print(ds)
        ti = ds.temp.values 
        ps = ds.ps.values 
        ti = np.squeeze(np.mean(ti,axis=0))
        ps = np.squeeze(np.mean(ps,axis=0))
        # ti  = np.squeeze(ti[-1,:,:,0])
        # ps = np.squeeze(ps[-1,:,0])

        t  = np.ones((50,64,128)) * ti[:,None,None]#[:,:,None]
        psurf = np.ones((64,128)) * ps[None,None]#[:,None]
        tsurf = np.ones((64,128)) * np.squeeze(ti[-1])[None,None]

        print(np.shape(tsurf), np.shape(np.squeeze(ti[-1])))

        u  = np.zeros_like(t)
        v  = np.zeros_like(t)
        sh = np.zeros_like(t)

        f = nc.Dataset(pss[j]+'bar_'+fracs[i]+'frac_ic.nc','w',format='NETCDF3_CLASSIC')

        f.createDimension('lon'  , 128) 
        f.createDimension('lat'  ,  64)
        f.createDimension('pfull',  50)

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

        longitude[:] = np.arange(0,360,360./128.)
        latitude[:]  = lat
        pfull[:]     = ds.pfull.values 
        temp[:,:,:]  = t
        ucomp[:,:,:] = u 
        vcomp[:,:,:] = v 
        sphum[:,:,:] = sh
        ps[:,:]      = psurf
        t_surf[:,:]  = tsurf

        f.close()