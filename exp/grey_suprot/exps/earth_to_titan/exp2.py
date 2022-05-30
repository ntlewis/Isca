"""
EXP2: 

res: T31
levels: high top 

omega: Earth 
a: Earth 
ps: Earth 

insolation: p2 
optical depth: frierson 


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
exp = Experiment('grey_T31HT_frierson_p2insol_om_1_1_a_1_1_ps_1bar_abs_0frac',
                  codebase=cb)


#### set constants and experiment arrays 


# horizontal resolution 
lons          = 96 
lats          = 48 
fourier       = 31
spherical     = 32

d_coeff       = 1. / (13.9 * 60. * 60.)

NCORES        = 24

timestep      = 300.

# vertical resolution 
level_type = 'high_top'





# planetary constants 
ps             = 100000.                       # surface pressure

grav           = 9.81                         # acceleration due to gravity 

omega          = 7.29e-5                      # rotation rate 

radius         = 6.371e6                      # radius (m)

R_d            = 287.                         # gas constant for dry air 

cp             = 1005.                        # specific heat capacity 
                                              # at constant pressure 

# convection constants 
dt_conv        = 3600. * 2.                   # timescale for dry convection 

gamma_conv     = 0.7                          # lapse rate for dry convection 

# ml code constants 
ml_depth       = 1.                           # mixed layer depth (m)

ml_albedo      = 0.                           # albedo for ml scheme 


# radiative transfer 
scheme         = 'frierson'                   # original grey rad or scatter?

insolation     = 'p2insol'                    # p2 insolation or perp eq?

solar_constant = 1360.                        # solar constant 

albedo         = 0.27                         # bond albedo 



tau_lw_s       = 2.                           # longwave optical depth (scat)

tau_sw_s       = 2./3. * 0.25 * 160. * tau_lw_s *12./16.  
                                              # shortwave optical depth 
                                              # in non-window region (scat)

pr             = 20000.                       # reference pressure for 
                                              # optical depths (unused)

lw_f           = 0.1                          # longwave p-dependence 
                                              # transition fraction (unused)

sw_window_frac = 1.                         # shortwave window fraction 
                                            # (if 1, atmosphere is transparent)

sw_window_flag = True                       # flag for shortwave window 





# set some things based on the choices above 

# vertical levels:
if level_type == 'low_top': 
    level_option        = 'hybrid'
    num_levels          = 25
    level_exponent      = 7.0 
    level_surf_res      = 0.5 
    level_scale_heights = 6.0
elif level_type == 'high_top': 
    level_option        = 'hybrid'
    num_levels          = 50 
    level_exponent      = 7.0 
    level_surf_res      = 0.5 
    level_scale_heights = 11.0




# radiative transfer: 
if scheme == 'frierson': 
    two_stream_gray = True
    two_stream_scatter = False 
elif scheme == 'generic': 
    two_stream_gray = False 
    two_stream_scatter = True 
if insolation == 'p2insol': 
    do_seasonal = False 
    perp_eq = False 
elif insolation == 'perp_eq': 
    if scheme == 'frierson':
        print('THIS WILL NOT WORK!!!!!')
        import sys 
        sys.exit()
    do_seasonal = False 
    perp_eq = True 















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
elif scheme == 'generic':
    diag.add_field('two_stream_scatter', 'olr', time_avg=True)
    diag.add_field('two_stream_scatter', 'swdn_toa', time_avg=True)
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
        'valid_range_t'            : [10.,800.],
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
        'two_stream_scatter': two_stream_scatter, 
        # use dry convection scheme 
        'convection_scheme' : 'dry', 
    },

  

    'two_stream_scatter_rad_nml': {
        'sw_optical_depth': 'generic', 
        'lw_abs_a'        : 0.0, 
        'lw_abs_b'        : 0.0, 
        'lw_abs_c'        : 0.0, 
        'lw_abs_d'        : 0.0,
        # longwave optical depth 
        'lw_abs_e'        : tau_lw_s, 
        'lw_abs_f'        : 0.0, 
        'lw_abs_pref'     : pr, #(unused)
        'lw_sca_a'        : 0.0,
        'lw_sca_b'        : 0.0,
        'sw_abs_a'        : 0.0,
        'sw_abs_b'        : 0.0, 
        # shortwave optical depth in non-window region 
        'sw_abs_c'        : tau_sw_s, 
        'sw_sca_a'        : 0.0, 
        'sw_sca_b'        : 0.0,
        # shortwave window fraction 
        'sw_window_frac'  : sw_window_frac, 
        'do_sw_window'    : sw_window_flag, 
        # assymetry factor for scattering (unused if sw_sca's are 0.0)
        'g_asym'          : 0.0, 
        # do true perpetual equinox, not frierson p2 
        'do_seasonal'     : do_seasonal,
        'perp_eq'         : perp_eq, 
        # apply albedo here, otherwise things get complicated with the window
        'solar_constant'  : solar_constant * (1 - albedo)
    },

        'two_stream_gray_rad_nml': {
        'rad_scheme'      : 'frierson',     
        'do_seasonal'     : do_seasonal,          
        'atm_abs'         : 0.0,      
        'solar_constant'  : solar_constant * (1 - albedo)
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
        'tconst':285., 
        # no evap because simulation is dry 
        'evaporation'  : False,   
        'depth'        : ml_depth, 
        # albedo zero because applied to solar constant instead                        
        'albedo_value' : ml_albedo,                          
    },

    'dry_convection_nml': {
        'tau'  : dt_conv,
        'Gamma': gamma_conv, 

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
        for i in range(2,11):
            exp.run(i, num_cores=NCORES, overwrite_data=False)
