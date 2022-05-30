""" 
This script configures a column model that uses Isca's columnwise physics routines. 

Useful for testing new convection / radiation parametrizations, as the dynamical core is 
bypassed so the model runs a gazillion times faster (especially if you're only simulating
one column). Can in principle simulate many (in lat and lon) at the same time. 

The wind is prescribed (it needs to be non-zero at the surface to allow for latent and 
sensible surface heat fluxes). Currently the user can set a namelist variable 'surface_wind'
that sets u_surf and v_surf = surface_wind / sqrt(2), so that wind_surf = sqrt(u_surf**2 + 
v_surf**2) = surface_wind. u and v at all other altitudes are set to zero (hardcoded). 

At the moment the model needs to use the vertical turbulent diffusion parameterization in order 
for the mixed layer code to work. This is not very consistent as the u and v wind are prescribed 
and so the u,v tendenency from the diffusion is thrown away. Hence an implicit assumption when 
using the column model is that 'the dynamics' would restore the surface winds to their prescribed 
speed, so that du/dt total is zero. 

The column model is currently initiated as a bit of a hack. The line 

'from isca import ColumnCodeBase'

sets a compiler flag -DCOLUMN_MODEL that tells the model to use the following files: 

atmos_column/column.F90
atmos_column/column_grid.F90
atmos_column/column_init_cond.F90
atmos_column/column_initialize_fields.F90

to initialize the model (including constructing the model grid), do the model timestepping 
(using a leapfrog scheme as before), and  handle input/output. 

Works with either hs_forcing, or the physics packages in idealized_moist_phys. Even when 
multiple columns are simulated, the model can only run on 1 core at the moment (will endeavour 
to fix this as some point). Also, the column model cannot read in topography input files. 

Any questions to Neil Lewis:  
neil.lewis@physics.ox.ac.uk 
"""


import os

import numpy as np

from isca import ColumnCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE


# column model only uses 1 core
NCORES = 1

# compile code 
base_dir = os.path.dirname(os.path.realpath(__file__))
cb = ColumnCodeBase.from_directory(GFDL_BASE)
cb.compile() 

# create an Experiment object to handle the configuration of model parameters
exp = Experiment('earth_grey_test_column', codebase=cb)

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('column', 'ps', time_avg=True)
diag.add_field('column', 'bk')
diag.add_field('column', 'pk')
diag.add_field('atmosphere', 'precipitation', time_avg=True)
diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('mixed_layer', 'flux_lhe', time_avg=True)
diag.add_field('mixed_layer', 'flux_t', time_avg=True)
#diag.add_field('column', 'sphum', time_avg=True)
diag.add_field('column', 'ucomp', time_avg=True)
diag.add_field('column', 'vcomp', time_avg=True)
diag.add_field('column', 'temp', time_avg=True)
diag.add_field('atmosphere', 'dt_ug_diffusion', time_avg=True)
diag.add_field('atmosphere', 'dt_vg_diffusion', time_avg=True)
diag.add_field('atmosphere', 'dt_qg_convection', time_avg=True)
diag.add_field('atmosphere', 'dt_tg_convection', time_avg=True)
diag.add_field('atmosphere', 'dt_tg_diffusion', time_avg=True)
diag.add_field('atmosphere', 'convection_rain', time_avg=True)
diag.add_field('two_stream_scatter', 'swdn_toa', time_avg=True)
diag.add_field('two_stream_scatter', 'swup_toa', time_avg=True)
diag.add_field('two_stream_scatter', 'olr', time_avg=True)
diag.add_field('two_stream_scatter', 'swdn_sfc', time_avg=True)
diag.add_field('two_stream_scatter', 'lwup_sfc', time_avg=True)
diag.add_field('two_stream_scatter', 'lwdn_sfc', time_avg=True)
diag.add_field('two_stream_scatter', 'net_lw_surf', time_avg=True)
diag.add_field('two_stream_scatter', 'tdt_rad', time_avg=True)
diag.add_field('two_stream_scatter', 'tdt_solar', time_avg=True)
diag.add_field('two_stream_scatter', 'flux_rad', time_avg=True)
diag.add_field('two_stream_scatter', 'flux_lw', time_avg=True)
diag.add_field('two_stream_scatter', 'flux_sw', time_avg=True)
diag.add_field('two_stream_scatter', 'flux_lw_up', time_avg=True)
diag.add_field('two_stream_scatter', 'flux_sw_up', time_avg=True)
diag.add_field('two_stream_scatter', 'flux_lw_down', time_avg=True)
diag.add_field('two_stream_scatter', 'flux_sw_down', time_avg=True)
diag.add_field('two_stream_scatter', 'flux_sw_down_direct', time_avg=True)
exp.diag_table = diag

#Empty the run directory ready to run
exp.clear_rundir()

#list planets parameters here 
ps = 100000
T_ini = 260.
grav = 9.81
tau_lw_s = 2. #was4 #(corresponds to 50,000 at surface) #473.56
tau_sw_s = 2./3. * 0.25 * 160. * tau_lw_s *12./16.  #.5#0.22 * tau_lw_s #0.22 is anti-greenhouse parameter 
solar_constant = 1360.
albedo = 0.27
omega = 7.29e-5
radius = 6.371e6
R_d = 287.
cp = 1005.
dt_conv = 3600. * 6. # slow timescale for dry convection 
ml_depth = 1. # meters 
pr = 20000. 
lw_f = 0.1 
sw_window_frac = 1. - 4./7.
sw_window_flag = True #True 






#Define values for the 'core' namelist
exp.namelist = namelist = Namelist({
    'main_nml':{
     'days'   : 2*3600,
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,
     'dt_atmos':2*14400,
     'current_date' : [1,1,1,0,0,0],
     'calendar' : 'thirty_day'
    },
    
    'atmosphere_nml': {
        'idealized_moist_model': True
    },

    'column_nml': {
        'lon_max': 1, # number of columns in longitude, default begins at lon=0.0
        'lat_max': 1, # number of columns in latitude, precise 
                      # latitude can be set in column_grid_nml if only 1 lat used. 
        'num_levels': 40,  # number of levels 
        'initial_sphum': 0.0, 
        'vert_coord_option': 'uneven_sigma',         # default: 'even_sigma'
        'scale_heights': 9.,#6.0,
        'exponent': 2.5,
        'surf_res': 0.1,
        'valid_range_t': [10., 3000.],         
        'reference_sea_level_press': ps # Pa 


    },

    'column_grid_nml': { 
        'lat_value': np.rad2deg(np.arcsin(1/np.sqrt(3))) # set latitude to that which causes insolation in frierson p2 radiation to be insolation / 4. 
        #'global_average': True # don't use this option at the moment
    },

    # set initial condition, NOTE: currently there is not an option to read in initial condition from a file. 
    'column_init_cond_nml': {
        'initial_temperature': T_ini, # initial atmospheric temperature 
        'surf_geopotential': 0.0, # applied to all columns 
        'surface_wind': 0. # as described above 
    },

    'idealized_moist_phys_nml': {
        'do_damping': False, # no damping in column model, surface wind prescribed 
        'turb':True,        # DONT WANT TO USE THIS, BUT NOT DOING SO IS STOPPING MIXED LAYER FROM WORKING
        'mixed_layer_bc':True, # need surface, how is this trying to modify the wind field? ****
        'do_simple': True, # simple RH calculation 
        'roughness_mom': 3.21e-05, # DONT WANT TO USE THIS, BUT NOT DOING SO IS STOPPING MIXED LAYER FROM WORKING
        'roughness_heat':3.21e-05,
        'roughness_moist':3.21e-05,                
        'two_stream_gray': False,     #Use grey radiation
        'two_stream_scatter':True, 
        'do_rrtm_radiation': False, 
        'do_socrates_radiation': False, 
        'convection_scheme': 'dry', #Use the simple Betts Miller convection scheme 
    },

    'two_stream_scatter_rad_nml': {
        'sw_optical_depth': 'generic', 
        'lw_abs_a': 0.0, 
        'lw_abs_b': 0.0, 
        'lw_abs_c': 0.0, 
        'lw_abs_d': 0.0,#tau_lw_s, #grav * tau_lw_s / ps, 
        'lw_abs_e': tau_lw_s, 
        'lw_abs_f': lw_f, 
        'lw_abs_pref': pr, 
        'lw_sca_a':0.0,
        'lw_sca_b':0.0,
        'sw_abs_a': 0.0,#grav * tau_sw_s / ps, 
        'sw_abs_b':0.0, 
        'sw_abs_c': tau_sw_s, 
        'sw_sca_a':0.0, 
        'sw_sca_b':0.0,
        'sw_window_frac':sw_window_frac, 
        'do_sw_window':sw_window_flag, 
        'g_asym' : 0.0, 
        'do_seasonal':False,
        'solar_constant': solar_constant * (1 - albedo) #2613.9 * (1 - 0.76)
    },

    'constants_nml': { # set constants to be those of Venus 
        'omega': omega, 
        'grav': grav, 
        'radius': radius, 
        'rdgas': R_d, 
        'PSTD': ps * 10., #dynes cm-2
        'PSTD_MKS': ps, #pa 
        'kappa': R_d / cp,
        #'CP_AIR': 900.,
    },


   # 'llcs_nml': { 
   #     'llcs_rhcrit': 0.8, 
   #     'llcs_timescale': dt_conv, 
   #     'cloud_option': 'ALL_RAIN', 
   #     'do_simple': True, 
   # },

    'dry_convection_nml': {
         'tau':dt_conv,
         'Gamma':0.65, 

     },
    
    'lscale_cond_nml': {
        'do_simple':True, # only rain 
        'do_evap':False,  # no re-evaporation of falling precipitation 
    },

    'surface_flux_nml': {
        'use_virtual_temp': False, # use virtual temperature for BL stability 
        'do_simple': True,
        'old_dtaudv': True    
    },

    'vert_turb_driver_nml': { # DONT WANT TO USE THIS, BUT NOT DOING SO IS STOPPING MIXED LAYER FROM WORKING
        'do_mellor_yamada': False,     # default: True
        'do_diffusivity': True,        # default: False
        'do_simple': True,             # default: False
        'constant_gust': 0.1,          # default: 1.0
        'use_tau': False
    },

    #Use a large mixed-layer depth, and the Albedo of the CTRL case in Jucker & Gerber, 2017
    'mixed_layer_nml': {
        'tconst' : T_ini,
        'prescribe_initial_dist':False,
        'evaporation':False,   
        'depth': ml_depth,                          #Depth of mixed layer used
        'albedo_value': 0.0,                  #Albedo value used             
    },
    
    'sat_vapor_pres_nml': {
        'do_simple':True, 
        'tcmin_simple' : -273, 
        'tcmax_simple' : 2000
    },

    # # define pressure coordinate 
    # 'vert_coordinate_nml': {
    #     'bk': [0.000000, 0.0117665, 0.0196679, 0.0315244, 0.0485411, 0.0719344, 0.1027829, 0.1418581, 0.1894648, 0.2453219, 0.3085103, 0.3775033, 0.4502789, 0.5244989, 0.5977253, 0.6676441, 0.7322627, 0.7900587, 0.8400683, 0.8819111, 0.9157609, 0.9422770, 0.9625127, 0.9778177, 0.9897489, 1.0000000],
    #     'pk': [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000],
    #    },

    # FMS Framework configuration
    'diag_manager_nml': {
        'mix_snapshot_average_fields': False  # time avg fields are labelled with time in middle of window
    },
    'fms_nml': {
        'domains_stack_size': 600000                        # default: 0
    },
    'fms_io_nml': {
        'threading_write': 'single',                         # default: multi
        'fileset_write': 'single',                           # default: multi
    },
})

#Lets do a run!
sw_window_frac_array = 1. - np.array([0.0,1./4.,2./4.,3./4.,1.]) 
names_array = ['0','1_4','1_2','3_4','1']

if __name__=="__main__":

    for n, sw_window_frac in enumerate(sw_window_frac_array):
        run_exp = exp.derive('earth_grey_test_column_frac_'+names_array[n])
        run_exp.namelist['two_stream_scatter_rad_nml']['sw_window_frac'] = sw_window_frac
        run_exp.run(1, use_restart=False, num_cores=NCORES, overwrite_data=True)
        #for i in range(2,16):
        #    run_exp.run(i, num_cores=NCORES, overwrite_data=False)
