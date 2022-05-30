"""
base_exp: 

res: T31
levels: high top 

omega: Earth / 16. (Titan)
a: Earth * 0.4 (Titan)
ps: 1.bar (Earth)
g: 9.81 ms-2 (Earth)

insolation: p2 (Titan params)
optical depth: no abs, lw robinson and catling  


"""

import os

import numpy as np

from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

###### SETTING UP THE MODEL #######
# set up base experiment including diagnostics and namelist 


base_dir = os.path.dirname(os.path.realpath(__file__))
cb = IscaCodeBase.from_directory(GFDL_BASE)
cb.compile()  

# template file name and experiment 
# template: 'grey_T31LT_frierson_p2insol_om_X_Y_a_X_Y_ps_Zbar_abs_X_Yfrac'
exp = Experiment('grey_T31HT_p2titan_om_1_16_a_4_10_ps_1bar_noabs',
                  codebase=cb)


#### set constants and experiment arrays 


# horizontal resolution 
lons          = 96 
lats          = 48 
fourier       = 31
spherical     = 32

d_coeff       = 1. / (9.42 * 60. * 60.)

NCORES        = 24

timestep      = 120.

# vertical resolution 
level_type = 'high_top'

t_ini = 300.



# planetary constants 
ps             = 100000.                      # surface pressure

grav           = 9.81                         # acceleration due to gravity 

omega          = 7.29e-5 / 16.                # rotation rate 

radius         = 6.371e6 * 0.4                # radius (m)

R_d            = 287.                         # gas constant for dry air 

cp             = 1005.                        # specific heat capacity 
                                              # at constant pressure 

# convection constants 
dt_conv    = 3600.
gamma_conv = 0.77


# ml code constants 
ml_depth       = 1.                           # mixed layer depth (m)

ml_albedo      = 0.                           # albedo for ml scheme 


# radiative transfer 
scheme         = 'titan'                   # original grey rad or scatter?

insolation     = 'p2insol'                    # p2 insolation or perp eq?



# robinson set 
tau_zero                = 4.324
pres_zero               = ps
do_rc_optical_depth     = True
tau_rc                  = (2.882 / 4.324) * tau_zero
pres_rc                 = (122453 / 150000) * ps
linear_tau              = 0.25
nonlin_tau_exponent     = 2.
sw_k1                   = 0.0 #35.168 
sw_k2                   = 0.0 #0.0599 

# these don't change 
bond_albedo             = 0.27 
sw_k1_frac              = 4./7. 
do_sw_split             = True 
solar_constant          = 1360.
del_sol                 = 0.88
del_sw                  = 0.0 
do_seasonal             = False 
solday                  = -10 
equinox_day             = 0 
dt_rad_avg              = -1 
use_time_average_coszen = False 
diabatic_acce           = 1.





# set some things based on the choices above 

# vertical levels:
if level_type == 'low_top': 
    level_option        = 'hybrid'
    num_levels          = 36
    level_exponent      = 7.0 
    level_surf_res      = 0.5 
    level_scale_heights = 5.0
elif level_type == 'high_top': 
    level_option        = 'hybrid'
    num_levels          = 50 
    level_exponent      = 7.0 
    level_surf_res      = 0.5 
    level_scale_heights = 11.0




# radiative transfer: 
if scheme == 'frierson': 
    two_stream_gray = True 
    two_stream_titan = False
elif scheme == 'titan': 
    two_stream_gray = False 
    two_stream_titan = True
if insolation == 'p2insol': 
    do_seasonal = False 
















# format for diagnostics 
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')

# which diagnostics to write (will get more out / higher frequency output for experiments of interest using future runs)
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'omega', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'height', time_avg=True)
diag.add_field('mixed_layer', 't_surf', time_avg=True)
if scheme == 'frierson':
    diag.add_field('two_stream', 'olr', time_avg=True)
    diag.add_field('two_stream', 'swdn_toa', time_avg=True)
elif scheme == 'titan':
    diag.add_field('titan_gray_rad', 'olr', time_avg=True)
    diag.add_field('titan_gray_rad', 'swdn_toa', time_avg=True)
    diag.add_field('titan_gray_rad', 'flux_lw', time_avg=True)
    diag.add_field('titan_gray_rad', 'flux_sw', time_avg=True)
#diag.add_field('damping', 'udt_rdamp', time_avg=True)

exp.diag_table = diag

#Empty the run directory ready to run
exp.clear_rundir()

#Define values for the 'core' namelist
exp.namelist = Namelist({
    'main_nml':{
     'days'   : 360,
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,
     'dt_atmos':timestep,
     'current_date' : [1,1,1,0,0,0],
     'calendar' : 'thirty_day'
    },

    'atmosphere_nml': {
        'idealized_moist_model': True
    },

    'spectral_dynamics_nml': {
        'damping_order'            : 2,                 
        'damping_coeff'            : d_coeff, 
        # laplacian sponge 
        'eddy_sponge_coeff'        : 1. / (0.4 * 60. * 60.),
        'water_correction_limit'   : 200.e2,
        'reference_sea_level_press': ps,
        'valid_range_t'            : [2.,2000.],
        'initial_sphum'            : [0.0],
        'robert_coeff'             : 0.03,
        'lon_max'                  : lons,
        'lat_max'                  : lats,
        'num_fourier'              : fourier,
        'num_spherical'            : spherical,
        'vert_coord_option'        : level_option,        
        'num_levels'               : num_levels,
        'exponent'                 : level_exponent, 
        'surf_res'                 : level_surf_res, 
        'scale_heights'            : level_scale_heights, 
    },



    'spectral_init_cond_nml': {
        'initial_temperature': t_ini, # initial atmospheric temperature 
    },



    'idealized_moist_phys_nml': {
        # no damping: sponge is laplacian in dynamical core 
        'do_damping'        : False, 
        'turb'              : True,
        'mixed_layer_bc'    : True,
        'do_virtual'        : False,
        'do_simple'         : True,
        'roughness_mom'     : 5.e-3,
        'roughness_heat'    : 5.e-3,
        'roughness_moist'   : 1.e-5,                
        'two_stream_gray'   : two_stream_gray,   
        'titan_gray'        : two_stream_titan, 
        # use dry convection scheme 
        'convection_scheme' : 'LLCS', 
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
        'diabatic_acce':diabatic_acce
    },

        'two_stream_gray_rad_nml': {
        'rad_scheme'      : 'frierson',     
        'do_seasonal'     : do_seasonal,          
        'atm_abs'         : 0.0,      
        'solar_constant'  : solar_constant * (1 - bond_albedo)
    },

    'vert_turb_driver_nml': {
        'do_mellor_yamada': False,    
        'do_diffusivity'  : True,       
        'do_simple'       : True,          
        'constant_gust'   : 1.0,         
        'use_tau'         : False
    },
    
    'diffusivity_nml': {
        'do_entrain': False,
        'do_simple' : True,
    },

    'surface_flux_nml': {
        'use_virtual_temp': False,
        'do_simple'       : True,
        'old_dtaudv'      : True    
    },

    'constants_nml':{
        # set planet constants 
        'omega'   : omega, # modified in experiment loop 
        'grav'    : grav, 
        'radius'  : radius, # modified in experiment loop 
        'rdgas'   : R_d, 
        'PSTD'    : ps * 10., 
        'PSTD_MKS': ps, 
        'kappa'   : R_d / cp,
    },



    'mixed_layer_nml': {
        # use initial condition (file provided in experiment loop)
        #'ic_neil_style': True, 
        'prescribe_initial_dist':True, 
        'tconst':t_ini, 
        # no evap because simulation is dry 
        'evaporation'  : False,   
        'depth'        : ml_depth, 
        # albedo zero because applied to solar constant instead                        
        'albedo_value' : ml_albedo,                          
    },

    'llcs_nml': {
        'alpha':gamma_conv,
        'llcs_timescale':dt_conv, 
        'do_simple':True
    },
    
    
    'lscale_cond_nml': {
        'do_simple': True,
        # no evaporation because experiments are dry 
        'do_evap': False 
    },
    
    'sat_vapor_pres_nml': {
        'do_simple'    : True,
        # range here doesn't matter as experiments are dry, but this prevents the model from crashing 
        'tcmin_simple' : -273, 
        'tcmax_simple' : 2000
    },




    # FMS Framework configuration
    'diag_manager_nml': {
        'mix_snapshot_average_fields': False  
    },

    'fms_nml': {
        'domains_stack_size': 600000                    
    },

    'fms_io_nml': {
        'threading_write': 'single',                       
        'fileset_write'  : 'single',                  
    },


})

# Run experiments 
if __name__=="__main__":


        exp.run(1, use_restart=False, num_cores=NCORES, overwrite_data=False)
        for i in range(2,81):
            exp.run(i, num_cores=NCORES, overwrite_data=False)
