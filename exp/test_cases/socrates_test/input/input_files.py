"""
Script that creates a file for prescribed polar heating based on Orlanski and Solman (2010) and Ruth Geen's H-S work
"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

def ideal_topo(y_min=25., y_max=65, h_0=4000., m=2., save_output=True):
    """
    #Bottom topography in the NH specified by setting the surface geopotential height.
    #Based on Sheshadri et al. (2015).
    Code based on Ruth Geen's hs_input_file.py script.
    """
    
    # Parameters
    # 1. y_min and y_max - centers topography at 45N
    # 2. m -  wave number of topography (1 or 2)
    # 3. h_0 = height of topography

    file = '/scratch/ntl203/Isca/input/dimensions_2d.nc'
    ds = xr.open_dataset(file, decode_times=False)
    
    template = ds.data * 0. + 1.
    topo_lat = (np.sin(np.pi * (ds.lat - y_min)/(y_max - y_min)))**2 * template
    topo_lon = np.cos(m * np.deg2rad(ds.lon)) * template

    zsurf = h_0 * topo_lat * topo_lon
    zsurf = zsurf.where(zsurf['lat']>y_min, 0)
    zsurf = zsurf.where(zsurf['lat']<y_max, 0)
    land_mask = template.where(zsurf != 0., 0)


    coord_list = ["lat", "lon"]
    zsurf = xr.Dataset(
         data_vars=dict(
             zsurf = (coord_list, zsurf.transpose('lat','lon').values),
             land_mask = (coord_list, land_mask.transpose('lat','lon').values)
         ),
         coords=ds.coords
    )
    
    if save_output:
        # NB filename should be 32 characters or less
        filename = 'h' + str(int(h_0)) + 'm' + str(int(m)) + 'l' + str(int(y_min)) + 'u' + str(int(y_max))
        print(len(filename))

        zsurf.to_netcdf(filename + '.nc', format="NETCDF3_CLASSIC",
             encoding = {'zsurf': {"dtype": 'float32', '_FillValue': None},
                    'land_mask': {"dtype": 'float32', '_FillValue': None},
                    "lat": {'_FillValue': None}, "lon": {'_FillValue': None},
                    "latb": {'_FillValue': None}, "lonb": {'_FillValue': None}}
                )
    
    return filename

def polar_heating(y_wid=15., th_mag=2., p_top = 500., p_th = 50., p_ref=800., save_output=True):
    
    # Parameter sweep
    # 1. Vary p_top - depth of forcing: 0:200:800 (1000 would be no forcing!)
    # 2. Vary th_mag - magnitude of forcing in K/s when centred on 800hPa  (0.5,1.,1.5,2.) Could try negative values too
    # 3. Vary y_wid - decay of forcing away from pole (10., 15., 20.)
    # 4. Vary p_th - sets vertical gradient of forcing at cap (25.,50.,75.) Sensitivity check that steepness of transition doesn't cause unexpected behaviour
    
    path = '/scratch/ntl203/Isca/input/'
    ozone_file = path+'rrtm_input_files/ozone_1990_notime.nc'
    data = xr.open_dataset(ozone_file, decode_times=False)
    
    template = data.ozone_1990 * 0. + 1.
    
    # Vary with latitude in similar way to Orlanski and Solman 2010
    #heat_lat = np.exp(-((data.lat - 90.)/y_wid)**2.) + np.exp(-((data.lat + 90.)/y_wid)**2.) * template
    heat_lat = np.exp(-((data.lat - 90.)/y_wid)**2.) * template # for north pole only
        
    # fix so that the function has magnitude 1 when p_top = p_ref, and otherwise scales to give constant net energy input
    # p_ref can be varied to alter the total input, which I think will be (1000-p_ref) * cp/g * th_mag
    
    if p_top==0.:  # If the heating is going right to the model top then make heating uniform in height, scaled to fit p_ref level
        polar_heat= th_mag * (1000. - p_ref)/1000. * heat_lat /86400.
    else:
        polar_heat = 0.5 * th_mag * (1000. - p_ref)/(1000. - p_top) * heat_lat * (1. + np.tanh((data.pfull - p_top)/p_th)) /86400.

    if int(y_wid) != 15:
        # scales strength of heating for different widths of heating vs. y_wid = 15 case
        # requires there to already be a file with all the same parameters except y_wid = 15
        print("scaling")
        file0 = 'w15' + 'a' + str(int(th_mag)) + 'p' + str(int(p_top)) + 'f' + str(int(p_ref)) + 'g' + str(int(p_th))
        ds0 = xr.open_dataset(path + 'polar_heating/' + file0 + '.nc')
        heat0 = heat1 = 0
        for i in range(len(data.lat)):
            heat0 += ds0.variables[file0][-1,i,0].data
            heat1 += polar_heat.transpose('pfull','lat','lon')[-1,i,0].data
        polar_heat = polar_heat / (heat1/heat0)

    coord_list = ["pfull", "lat", "lon"]
    polar_heat = xr.Dataset(
         data_vars=dict(
             polar_heat = (coord_list, polar_heat.transpose('pfull','lat','lon').values)
         ),
         coords=data.coords
    )

    if save_output:
        # NB filename should be 32 characters or less
        filename = 'w' + str(int(y_wid)) + 'a' + str(int(th_mag)) + 'p' + str(int(p_top)) + 'f' + str(int(p_ref)) + 'g' + str(int(p_th))
        print(len(filename))
        polar_heat = polar_heat.rename({"polar_heat" : filename})
        
        polar_heat.to_netcdf(filename + '.nc', format="NETCDF3_CLASSIC",
             encoding = {filename: {"dtype": 'float32', '_FillValue': None},
                    "lat": {'_FillValue': None}, "lon": {'_FillValue': None},
                    "latb": {'_FillValue': None}, "lonb": {'_FillValue': None},
                    "pfull": {'_FillValue': None}, "phalf": {'_FillValue': None}}
                )    

    return filename

def heat_perturb(q_0=6, m=2, y_cen=45, p_0=800, p_t=200, save_output=True):
    """
    Tropospheric diabatic heating perturbation to induce NH winter-like wave activity.
    Based on Lindgren et al. (2018).
    Code based on Ruth Geen's hs_input_file.py script.
    """
    # Parameters
    # 1. q_0 - magnitude of forcing in K/day
    # 2. m - longitudinal wave number (1 or 2)
    # 3. y_cen - center of heating in latitude
    # 4. p_0 and p_t - lower and upper bounds in pressure (hPa) respectively
    y_wid = 0.175*360/(2*np.pi)

    file = '/scratch/ntl203/Isca/input/dimensions_3d.nc'
    ds = xr.open_dataset(file, decode_times=False)
    
    template = ds.data * 0. + 1.
    heat_lat = np.exp(-0.5 * ((ds.lat - y_cen)/y_wid)**2.) * template
    heat_lon = np.cos(m * np.deg2rad(ds.lon)) * template
    heat_p = np.sin(np.pi * np.log(ds.pfull/p_0)/np.log(p_t/p_0)) * template

    heat = q_0 * heat_lat * heat_lon *  heat_p /86400.  # convert to K/s
    heat = heat.where(heat['pfull']<p_0, 0)
    heat = heat.where(heat['pfull']>p_t, 0)

    coord_list = ["pfull", "lat", "lon"]
    heat = xr.Dataset(
         data_vars=dict(
             heat = (coord_list, heat.transpose('pfull','lat','lon').values)
         ),
         coords=ds.coords
    )
    
    if save_output:
        # NB filename should be 32 characters or less
        filename = 'q' + str(int(q_0)) + 'm' + str(int(m)) + 'y' + str(int(y_cen)) + 'l' + str(int(p_0)) + 'u' + str(int(p_t))
        print(len(filename))
        heat = heat.rename({"heat" : filename})

        heat.to_netcdf(filename + '.nc', format="NETCDF3_CLASSIC",
             encoding = {filename: {"dtype": 'float32', '_FillValue': None},
                    "lat": {'_FillValue': None}, "lon": {'_FillValue': None},
                    "latb": {'_FillValue': None}, "lonb": {'_FillValue': None},
                    "pfull": {'_FillValue': None}, "phalf": {'_FillValue': None}}
                )
    
    return filename

def combo_heat1(y_wid=15., th_mag=4., p_top = 800., p_th = 50., p_ref=800., q_0=6., m=2., y_cen=45., p_0=800., p_t=200., save_output=True):
    
    path = '/scratch/ntl203/Isca/input/'
    file = path+'dimensions_3d.nc'
    ds = xr.open_dataset(file, decode_times=False)
    template = ds.data * 0. + 1.

    # Start with polar heating
    # Parameters
    # 1. p_top - depth of forcing: 0:200:800 (1000 would be no forcing!)
    # 2. th_mag - magnitude of forcing in K/s when centred on 800hPa  (0.5,1.,1.5,2.) Could try negative values too
    # 3. y_wid - decay of forcing away from pole (10., 15., 20.)
    # 4. p_th - sets vertical gradient of forcing at cap (25.,50.,75.) Sensitivity check that steepness of transition doesn't cause unexpected behaviour
    
    # Vary with latitude in similar way to Orlanski and Solman 2010
    #heat_lat = np.exp(-((data.lat - 90.)/y_wid)**2.) + np.exp(-((data.lat + 90.)/y_wid)**2.) * template
    heat_lat = np.exp(-((ds.lat - 90.)/y_wid)**2.) * template # for north pole only
        
    # fix so that the function has magnitude 1 when p_top = p_ref, and otherwise scales to give constant net energy input
    # p_ref can be varied to alter the total input, which I think will be (1000-p_ref) * cp/g * th_mag
    
    if p_top==0.:  # If the heating is going right to the model top then make heating uniform in height, scaled to fit p_ref level
        polar_heat = th_mag * (1000. - p_ref)/1000. * heat_lat /86400.
    else:
        polar_heat = 0.5 * th_mag * (1000. - p_ref)/(1000. - p_top) * heat_lat * (1. + np.tanh((ds.pfull - p_top)/p_th)) /86400.
        
    if int(y_wid) != 15:
        # scales strength of heating for different widths of heating vs. y_wid = 15 case
        # requires there to already be a file with all the same parameters except y_wid = 15
        print("scaling")
        file0 = 'w15' + 'a' + str(int(th_mag)) + 'p' + str(int(p_top)) + 'f' + str(int(p_ref)) + 'g' + str(int(p_th))
        ds0 = xr.open_dataset(path + 'polar_heating/' + file0 + '.nc')
        heat0 = heat1 = 0
        for i in range(len(ds.lat)):
            heat0 += ds0.variables[file0][-1,i,0].data
            heat1 += polar_heat.transpose('pfull','lat','lon')[-1,i,0].data
        polar_heat = polar_heat / (heat1/heat0)

    # Now do zonally asymmetric heat perturbation
    # Parameters
    # 1. q_0 - magnitude of forcing in K/day
    # 2. m - longitudinal wave number (1 or 2)
    # 3. y_cen - center of heating in latitude
    # 4. p_0 and p_t - lower and upper bounds in pressure (hPa) respectively
    y_w = 0.175*360/(2*np.pi)

    heat_lat = np.exp(-0.5 * ((ds.lat - y_cen)/y_w)**2.) * template
    heat_lon = np.cos(m * np.deg2rad(ds.lon)) * template
    heat_p = np.sin(np.pi * np.log(ds.pfull/p_0)/np.log(p_t/p_0)) * template

    heat = q_0 * heat_lat * heat_lon *  heat_p /86400.  # convert to K/s
    heat = heat.where(heat['pfull']<p_0, 0)
    heat = heat.where(heat['pfull']>p_t, 0)

    combo_heat = polar_heat + heat

    coord_list = ["pfull", "lat", "lon"]
    combo_heat = xr.Dataset(
         data_vars=dict(
             combo_heat = (coord_list, combo_heat.transpose('pfull','lat','lon').values)
         ),
         coords=ds.coords
    )

    if save_output:
        # NB filename should be 32 characters or less
        filename = 'w' + str(int(y_wid)) + 'a' + str(int(th_mag)) + 'p' + str(int(p_top)) + 'f' + str(int(p_ref)) + 'g' + str(int(p_th)) +\
            '_' 'q' + str(int(q_0)) + 'm' + str(int(m)) + 'y' + str(int(y_cen)) + 'l' + str(int(p_0)) + 'u' + str(int(p_t))
        print(len(filename))
        combo_heat = combo_heat.rename({"combo_heat" : filename})
        
        combo_heat.to_netcdf(filename + '.nc', format="NETCDF3_CLASSIC",
             encoding = {filename: {"dtype": 'float32', '_FillValue': None},
                    "lat": {'_FillValue': None}, "lon": {'_FillValue': None},
                    "latb": {'_FillValue': None}, "lonb": {'_FillValue': None},
                    "pfull": {'_FillValue': None}, "phalf": {'_FillValue': None}}
                )    

    return filename

def offpole_heating(q_0=4., x_cen=75., y_cen=180., x_wid=5., y_wid=30., p_top = 800., p_th = 50., p_ref=800., save_output=True):
    
    # Parameters
    # 1. q_0 - magnitude of forcing in K/day
    # 2. x/y_cen - center of heating in latitude/longitude
    # 3. x/y_wid - width of heating in latitude/longitude
    # 4. Vary p_top - depth of forcing
    # 4. Vary p_th - sets vertical gradient of forcing
    
    path = '/scratch/ntl203/Isca/input/'
    file = path+'dimensions_3d.nc'
    ds = xr.open_dataset(file, decode_times=False)
    template = ds.data * 0. + 1.
    
    heat_lat = np.exp(-0.5 * ((ds.lat - x_cen)/x_wid)**2.) * template
    heat_lon = np.exp(-0.5 * ((ds.lon - y_cen)/y_wid)**2.) * template
    heat_p = 0.5 * (1. + np.tanh((ds.pfull - p_top)/p_th)) * template

    heat = q_0 * (1000. - p_ref)/(1000. - p_top) * heat_lat * heat_lon *  heat_p /86400.  # convert to K/s

    coord_list = ["pfull", "lat", "lon"]
    heat = xr.Dataset(
         data_vars=dict(
             heat = (coord_list, heat.transpose('pfull','lat','lon').values)
         ),
         coords=ds.coords
    )

    if save_output:
        # NB filename should be 32 characters or less
        filename = 'a' + str(int(q_0)) + 'x' + str(int(x_cen)) + 'y' + str(int(y_cen)) +\
            'w' + str(int(x_wid)) + 'v' + str(int(y_wid)) + 'p' + str(int(p_top))
        print(len(filename))
        heat = heat.rename({"heat" : filename})
        
        heat.to_netcdf(filename + '.nc', format="NETCDF3_CLASSIC",
             encoding = {filename: {"dtype": 'float32', '_FillValue': None},
                    "lat": {'_FillValue': None}, "lon": {'_FillValue': None},
                    "latb": {'_FillValue': None}, "lonb": {'_FillValue': None},
                    "pfull": {'_FillValue': None}, "phalf": {'_FillValue': None}}
                )

    return filename

def plot_vertical(folder, filename):

    file = filename + '.nc'
    ds = xr.open_dataset(file)
    if folder == 'polar_heating':
        h = int(filename.partition("a")[2][0])/86400
    elif folder == 'asymmetry':
        h = int(filename.partition("q")[2][0])/86400
    inc = 0.5e-5

    lat = ds.coords['lat'].data
    p = ds.coords['pfull'].data
    heat = ds.sel(lon=180, method='nearest').variables[filename]
    
    # Plot
    plt.figure(figsize=(8,6))
    cs = plt.contourf(lat, p, heat, cmap='Reds', levels=np.arange(0, h+inc, inc))
    cb = plt.colorbar(cs)
    cb.set_label(label=r'Heating (K s$^{-1}$)', size='x-large')
    cb.ax.tick_params(labelsize='x-large')
    plt.xlabel(r'Latitude ($\degree$N)', fontsize='x-large')
    plt.xticks([30, 50, 70, 90], ['30', '50', '70', '90'])
    plt.xlim(20, 90)
    plt.ylim(max(p), 100)
    plt.yscale('log')
    plt.ylabel('Pressure (hPa)', fontsize='x-large')
    #plt.title(filename, fontsize='x-large')
    plt.tick_params(axis='both', labelsize = 'x-large', which='both', direction='in')
    plt.savefig(filename+'v.pdf', bbox_inches = 'tight')

    return plt.close()

def plot_horizontal1(folder, filename):
    
    file = filename + '.nc'
    ds = xr.open_dataset(file)
    if folder == 'polar_heating':
        h = int(filename.partition("a")[2][0])/86400
    elif folder == 'asymmetry':
        h = int(filename.partition("q")[2][0])/86400
    inc = 0.5e-5

    lat = ds.coords['lat'].data
    lon = ds.coords['lon'].data
    heat = ds.sel(pfull=1000, method='nearest').variables[filename]

    #Plot
    plt.figure(figsize=(8,6))
    plt.contourf(lon, lat, heat, cmap='RdBu_r', levels=np.arange(0, h+inc, inc))
    plt.colorbar(label=r'Heating (K s$^{-1}$)')
    plt.xlabel(r'Longitude ($\degree$)', fontsize='x-large')
    plt.ylabel(r'Latitude ($\degree$)', fontsize='x-large')
    plt.title(filename, fontsize='x-large')
    plt.tick_params(axis='both', labelsize = 'x-large', which='both', direction='in')
    plt.savefig(filename+'h1.pdf', bbox_inches = 'tight')

    return plt.close()

def plot_horizontal2(folder, filename):
    
    file = filename + '.nc'
    ds = xr.open_dataset(file)
    if folder == 'polar_heating':
        h = int(filename.partition("a")[2][0])/86400
    elif folder == 'asymmetry':
        h = int(filename.partition("q")[2][0])/86400
    inc = 0.25e-5

    lat = ds.coords['lat'].data
    lon = ds.coords['lon'].data
    heat = ds.sel(pfull=1000, method='nearest').variables[filename]

    #Plot
    plt.figure(figsize=(8,6))
    ax = plt.axes(projection=ccrs.NorthPolarStereo())
    cs = plt.contourf(lon, lat, heat, cmap='RdBu_r', levels=np.arange(0, h+inc, inc), transform = ccrs.PlateCarree())
    cb = plt.colorbar(cs, pad=0.1)
    cb.set_label(label=r'Heating (K s$^{-1}$)', size='x-large')
    cb.ax.tick_params(labelsize='x-large')
    ax.coastlines()
    ax.set_global()
    ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
    ax.set_extent([-180, 180, 0, 90], crs=ccrs.PlateCarree())
    plt.savefig(filename+'h2.pdf', bbox_inches = 'tight')

    return plt.close()

if __name__ == '__main__': 
    option = input("a) polar heat, b) heat perturb, c) idealised topography, d) combo of a and b?, e) off-pole heat, or f) combo of e and b?")

    H = 8
    p0 = 1000

    if option =='a':
        filename = polar_heating()
        plot_vertical('polar_heating', filename)
    elif option =='b':
        filename = heat_perturb()
        plot_horizontal('asymmetry', filename)
        plot_vertical('asymmetry', filename)
    elif option =='c':
        filename = ideal_topo()
        #plot_horizontal('asymmetry', filename)
    elif option =='d':
        filename = combo_heat1()
        plot_vertical('asymmetry', filename)
        plot_horizontal('asymmetry', filename)
    elif option =='e':
        filename = offpole_heating()
        plot_horizontal1('polar_heating', filename)
        plot_horizontal2('polar_heating', filename)
        plot_vertical('polar_heating', filename)
