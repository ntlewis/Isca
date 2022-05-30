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
exp = Experiment('venus_col_logcheck_finalFlux_52H145_leblevs_oldlapse', codebase=cb)

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 360, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('column', 'ps', time_avg=True)
diag.add_field('column', 'bk')
diag.add_field('column', 'pk')
diag.add_field('atmosphere', 'precipitation', time_avg=True)
diag.add_field('mixed_layer', 't_surf', time_avg=True)
#diag.add_field('mixed_layer', 'flux_lhe', time_avg=True)
diag.add_field('mixed_layer', 'flux_t', time_avg=True)
#diag.add_field('column', 'sphum', time_avg=True)
diag.add_field('column', 'ucomp', time_avg=True)
diag.add_field('column', 'vcomp', time_avg=True)
diag.add_field('column', 'temp', time_avg=True)
diag.add_field('column', 'height', time_avg=True)
diag.add_field('atmosphere', 'dt_ug_diffusion', time_avg=True)
diag.add_field('atmosphere', 'dt_vg_diffusion', time_avg=True)
# diag.add_field('atmosphere', 'dt_qg_convection', time_avg=True)
diag.add_field('atmosphere', 'dt_tg_convection', time_avg=True)
diag.add_field('atmosphere', 'dt_tg_diffusion', time_avg=True)
# diag.add_field('atmosphere', 'convection_rain', time_avg=True)
diag.add_field('titan_gray_rad', 'olr', time_avg=True)
diag.add_field('titan_gray_rad', 'swdn_sfc', time_avg=True)
diag.add_field('titan_gray_rad', 'swdn_toa', time_avg=True)
diag.add_field('titan_gray_rad', 'net_lw_surf', time_avg=True)
diag.add_field('titan_gray_rad', 'lwdn_sfc', time_avg=True)
diag.add_field('titan_gray_rad', 'lwup_sfc', time_avg=True)
diag.add_field('titan_gray_rad', 'tdt_rad', time_avg=True)
diag.add_field('titan_gray_rad', 'tdt_solar', time_avg=True)
diag.add_field('titan_gray_rad', 'flux_rad', time_avg=True)
diag.add_field('titan_gray_rad', 'flux_lw', time_avg=True)
diag.add_field('titan_gray_rad', 'flux_sw', time_avg=True)
diag.add_field('titan_gray_rad', 'coszen', time_avg=True)
diag.add_field('titan_gray_rad', 'fracsun', time_avg=True)
diag.add_field('titan_gray_rad', 'tau', time_avg=True)
exp.diag_table = diag

#Empty the run directory ready to run
exp.clear_rundir()

#list planets parameters here 

# planetary parameters 
ps = 92.e4 #100000.#150000
grav = 8.87#1.35#9.81#1.35
omega = 2.99e-7#4.56e-6
radius = 6.051e6#2.576e6
R_d = 188.#290.# 287. #290.
cp = 900.#1000.#1005.

# init parameters 
T_ini = 710.#290.#100.

# convection parameters 
gamma_conv = .77#98
dt_conv = 3600. #* 8
do_relax = False#True#False

# mixed layer parameters
ml_depth = 1. # meters 
ml_albedo = 0.

# rt parameters 
var_cp = False#True#False 
# # mckay set 
tau_zero                = 1250#1100#5.#4.54#4.47 
pres_zero               = ps##150000
do_rc_optical_depth     = False 
tau_rc                  = -1. 
pres_rc                 = -1.
linear_tau              = 0.0
nonlin_tau_exponent     = 1.#4./3.
lin_tau_exponent        = 1.#2./3.

sw_k1                   = 0.0#141.22/4.#35.167 
sw_k2                   = 0.0#0.0643#0.0613 
do_log_sw = True
do_sw_venus = False
logtau_B  = 0.1074#0.11204#0.11069#0.09980#0.09229
logtau_sw_s =  0.899#0.90057#1.5248#1.497#1.897
logtau_0 = 1143.799#903.88#77.57#40.128#148.9#106#193.7#356.7
logtau_siglev = 0.004#2.47255e-3#0.00022#0.0005761
logtau_pres_zero = ps
logf = 1 - 0.7 #25
sw_tau_exponent = lin_tau_exponent

# robinson set 
# tau_zero                = 4.324 #2. #4.324
# pres_zero               = ps
# do_rc_optical_depth     = True
# tau_rc                  = (2.882 / 4.324) * tau_zero
# pres_rc                 = (122453 / 150000) * ps
# linear_tau              = 0.25
# nonlin_tau_exponent     = 2.
# sw_k1                   = 35.168 #0.0 #35.168 #0.0 #
# sw_k2                   = 0.0599 #0.0 #0.0599 # #

# these don't change 
bond_albedo             = 0.0#27 
sw_k1_frac              = 1.#4./7. #1.
do_sw_split             = False#True # False
solar_constant          = 158.16*4#3.6 * 4. #1360. #
del_sol                 = 1.4 
del_sw                  = 0.0 
do_seasonal             = False 
solday                  = -10 
equinox_day             = 0 
dt_rad_avg              = -1 
use_time_average_coszen = False 
diabatic_acce           = 1.



pk = np.flip(np.array([0.,656.,2.61e3, 7.67e3, 1.73e4, 3.3e4, 5.57e4, 8.63e4, 
                       1.25e5, 1.7e5, 2.21e5, 2.75e5, 3.27e5, 3.74e5, 4.13e5, 
                       4.4e5, 4.52e5, 4.5e5, 4.35e5, 4.01e5, 3.29e5, 2.37e5, 
                       1.7e5, 1.22e5, 8.73e4, 6.26e4, 4.48e4, 3.21e4, 2.3e4, 
                       1.65e4, 1.18e4, 8.47e3, 6.07e3, 4.35e3, 3.12e3, 2.23e3, 
                       1.6e3, 1.15e3, 8.21e2, 5.89e2, 4.22e2, 3.02e2, 2.12e2, 
                       1.55e2, 1.11e2, 79.7, 51.5, 28.3, 14.3, 6.11, 0.])) * ps / 92.e5 
bk = np.flip(np.concatenate([np.array([1., 9.99e-1, 9.95e-1, 9.85e-1, 9.66e-1, 
                               9.35e-1, 8.91e-1, 8.33e-1, 7.62e-1, 6.79e-1, 5.86e-1, 
                               4.89e-1, 3.91e-1, 2.99e-1, 2.15e-1, 1.43e-1, 8.75e-2, 
                               4.9e-2, 2.29e-2, 6.66e-3, 2.67e-4, 5.38e-8, 2.69e-15]), 
                                np.zeros(28)]))

new_pk = pk.tolist() #(pk + ps * bk).tolist()
new_bk = bk.tolist() #(np.zeros_like(bk)).tolist()



#Define values for the 'core' namelist
exp.namelist = namelist = Namelist({
    'main_nml':{
     'days'   : 3600,
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,
     'dt_atmos': 120,##1800,
     'current_date' : [1,1,1,0,0,0],
     'calendar' : 'thirty_day'
    },
    
    'atmosphere_nml': {
        'idealized_moist_model': True
    },

    # 'column_nml': {
    #     'lon_max': 1, # number of columns in longitude, default begins at lon=0.0
    #     'lat_max': 1, # number of columns in latitude, precise 
    #                   # latitude can be set in column_grid_nml if only 1 lat used. 
    #     'num_levels': 50,  # number of levels 
    #     'initial_sphum': 0.0, 
    #     'vert_coord_option': 'uneven_sigma',         # default: 'even_sigma'
    #     'scale_heights': 11.0,#6.0,
    #     'exponent': 7.0,
    #     'surf_res': 0.5 ,
    #     'valid_range_t': [10., 3000.],         
    #     'reference_sea_level_press': ps # Pa 


    # },

    'column_nml': {
        'lon_max': 1, # number of columns in longitude, default begins at lon=0.0
        'lat_max': 1, # number of columns in latitude, precise 
                      # latitude can be set in column_grid_nml if only 1 lat used. 
        'num_levels': 50,  # number of levels 
        'initial_sphum': 0.0, 
        'vert_coord_option': 'input',         # default: 'even_sigma'
#         'scale_heights': 14.5,#6.0,
#         'exponent': 3.0,
#         'surf_res': 0.15,
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
        'two_stream_scatter':False, 
        'titan_gray':True, 
        'do_rrtm_radiation': False, 
        'do_socrates_radiation': False, 
        'convection_scheme': 'DRYADJ', #Use the simple Betts Miller convection scheme 
    },

    'titan_gray_rad_nml': {
        'tau_rc': tau_rc, 
        'pres_rc': pres_rc, 
        'do_rc_optical_depth': do_rc_optical_depth,
        'tau_zero': tau_zero, 
        'pres_zero': pres_zero, 
        'linear_tau':linear_tau, 
        'nonlin_tau_exponent': nonlin_tau_exponent, 
        'bond_albedo': bond_albedo, 
        'sw_k1': sw_k1,
        'sw_k2': sw_k2, 
        'sw_k1_frac': sw_k1_frac, 
        'do_sw_split': do_sw_split,
        'solar_constant':solar_constant,
        'del_sol': del_sol,
        'del_sw':  del_sw, 
        'do_seasonal': do_seasonal, 
        'solday':solday, 
        'equinox_day':equinox_day,
        'dt_rad_avg':dt_rad_avg, 
        'use_time_average_coszen':use_time_average_coszen, 
        'diabatic_acce':diabatic_acce, 
        'do_log_sw3':do_log_sw, 
        'logtau_B': logtau_B, 
        'logtau_sw_s':logtau_sw_s, 
        'logtau_0': logtau_0, 
        'logtau_siglev':logtau_siglev,
        'logtau_pres_zero': logtau_pres_zero, 
        'logf':logf, 
        'lin_tau_exponent': lin_tau_exponent, 
        'sw_tau_exponent': sw_tau_exponent, 
        'var_cp' : var_cp, 
        'do_sw_venus':do_sw_venus
    },

    'constants_nml': { # set constants to be those of Venus 
        'omega': omega, 
        'grav': grav, 
        'radius': radius, 
        'rdgas': R_d, 
        'PSTD': ps * 10., #dynes cm-2
        'PSTD_MKS': ps, #pa 
        'kappa': R_d / cp,
    },




    'dry_adj_nml': {
        'alpha':gamma_conv,
        'tau_relax':dt_conv,#, 
        'do_relax':do_relax, 
        #'itermax':20,
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
        'do_diffusivity': False,        # default: False
        'do_simple': True,             # default: False
        'constant_gust': 0.0,          # default: 1.0
        'use_tau': False
    },

    #Use a large mixed-layer depth, and the Albedo of the CTRL case in Jucker & Gerber, 2017
    'mixed_layer_nml': {
        'tconst' : T_ini,
        'prescribe_initial_dist':False,
        'evaporation':False,   
        'depth': ml_depth,                          #Depth of mixed layer used
        'albedo_value': ml_albedo,                  #Albedo value used             
    },
    
    'sat_vapor_pres_nml': {
        'do_simple':True, 
        'tcmin_simple' : -273, 
        'tcmax_simple' : 2000
    },

    # define pressure coordinate 
    'vert_coordinate_nml': {
        'bk': new_bk,
        'pk': new_pk,
       },

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
if __name__=="__main__":
    exp.run(1, use_restart=False, num_cores=NCORES, overwrite_data=True)
    for i in range(2,101):
        exp.run(i, num_cores=NCORES, overwrite_data=True)
