"""
base_exp: 

res: T31
levels: high top 

omega: Earth / 16. (Titan)
a: Earth * 0.4 (Titan)
ps: 1.45bar (Titan)
g: 1.35 ms-2 (Titan)

insolation: p2 (Titan params)
optical depth: robinson and catling  


"""

import os

import numpy as np

from isca import GreyCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

###### SETTING UP THE MODEL #######
# set up base experiment including diagnostics and namelist 


base_dir = os.path.dirname(os.path.realpath(__file__))
cb = GreyCodeBase.from_directory(GFDL_BASE)
cb.compile()  

# template file name and experiment 
# template: 'grey_T31LT_frierson_p2insol_om_X_Y_a_X_Y_ps_Zbar_abs_X_Yfrac'
exp = Experiment('grey_T31L48_p2insol_om_1_1_a_1_8_ps_1_spg_gamma_7_10',
                  codebase=cb)


#### set constants and experiment arrays 


# horizontal resolution 
lons          = 96 
lats          = 48
fourier       = 31
spherical     = 32

d_coeff       = 1. / (5.8 * 60. * 60.) # 12.2

NCORES        = 24

timestep      = 300.

# vertical resolution 
#level_type = 'high_top'

t_ini = 300.



# planetary constants 
ps             = 100000.                      # surface pressure

grav           = 9.81                         # acceleration due to gravity 

omega          = 7.29e-5 #/ 16.               # rotation rate 

radius         = 6.371e6 * 0.125              # radius (m)

R_d            = 287.                         # gas constant for dry air 

cp             = 1005.                        # specific heat capacity 
                                              # at constant pressure 

# convection constants 
dt_conv    = 7200.
gamma_conv = 0.7


# ml code constants 
ml_depth       = 1.                           # mixed layer depth (m)

ml_albedo      = 0.                           # albedo for ml scheme 


# radiative transfer 
scheme         = 'titan'                   # original grey rad or scatter?

insolation     = 'p2insol'                    # p2 insolation or perp eq?



# robinson set 
tau_zero                = 1.15
pres_zero               = ps
do_rc_optical_depth     = True
tau_rc                  = (0.11/1.15)*tau_zero
pres_rc                 = (30403./1.e5)*pres_zero
linear_tau              = 0.25
nonlin_tau_exponent     = 2.
sw_k1                   = 0.0 #51.8 #35.168 #0.0 #
sw_k2                   = 0.0 #0.071 #0.0599 # #
sw_k3                   = 0.0 #8.0745  / 1.305e2 * sw_k1 
sw_gamma                = 3.395e-1

# these don't change 
bond_albedo             = 0.27 
sw_k1_frac              = 4./7. 
do_sw_split             = True 
solar_constant          = 1360.
del_sol                 = 1.2
del_sw                  = 0.0 
do_seasonal             = False 
solday                  = -10 
equinox_day             = 0 
dt_rad_avg              = -1 
use_time_average_coszen = False 
diabatic_acce           = 1.




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
        #'eddy_sponge_coeff'        : 1. / (0.4 * 60. * 60.),
        'water_correction_limit'   : 200.e2,
        'reference_sea_level_press': ps,
        'valid_range_t'            : [2.,2000.],
        'initial_sphum'            : [0.0],
        'robert_coeff'             : 0.03,
        'lon_max'                  : lons,
        'lat_max'                  : lats,
        'num_fourier'              : fourier,
        'num_spherical'            : spherical,
        'vert_coord_option'        : 'input',        
        'num_levels'               : 48,
    },



    'spectral_init_cond_nml': {
        'initial_temperature': t_ini, # initial atmospheric temperature 
    },



    'idealized_moist_phys_nml': {
        # no damping: sponge is laplacian in dynamical core 
        'do_damping'        : True, 
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

    'damping_driver_nml': {
        'do_rayleigh': True,
        'trayfric': -(1./2.),              # neg. value: time in *days*
        'sponge_pbottom':  1000.,           #Bottom of the model's sponge down to 50Pa (units are Pa)
        'do_conserve_energy': True,            
        'rayleigh_eddy_only':False, 
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

    'vert_coordinate_nml':{
        'bk': [0.00000,       0.00000,       0.00000,      
               0.00000,       0.00000,       0.00000,      
               0.00000,       0.00000,       0.00000,      
               0.00000,       0.00000,       0.00000,      
               0.00000,       0.00000,       0.00000,      
               0.00000,       0.00000,       0.00000,      
               0.00000,       0.00000,       0.00000,      
               0.00000,       0.00000,       0.00000,      
               0.00000,       0.00000,       0.00000,      
               0.00000,       0.00000,       0.00000,      
               0.00000,       0.00813,       0.03224,      
               0.07128,       0.12445,       0.19063,      
               0.26929,       0.35799,       0.45438,      
               0.55263,       0.64304,       0.71703,      
               0.77754,       0.82827,       0.87352,      
               0.91502,       0.95235,       0.98511,      
               1.00000  ], 
        'pk': [40.00000,     100.00000,     200.00000,     
              350.00000,     550.00000,     800.00000,     
             1085.00000,    1390.00000,    1720.00000,     
             2080.00000,    2470.00000,    2895.00000,     
             3365.00000,    3890.00000,    4475.00000,     
             5120.00000,    5830.00000,    6608.00000,     
             7461.00000,    8395.00000,    9424.46289,     
            10574.46900,   11864.80330,   13312.58850,     
            14937.03770,   16759.70760,   18804.78670,     
            21099.41250,   23674.03720,   26562.82650,     
            29804.11680,   32627.31601,   34245.89759,     
            34722.29104,   34155.20062,   32636.50533,     
            30241.08406,   27101.45052,   23362.20912,     
            19317.04955,   15446.17194,   12197.45091,     
             9496.39912,    7205.66920,    5144.64339,     
             3240.79521,    1518.62245,       0.00000,     
                0.00000]
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
        for i in range(2,6):
            exp.run(i, num_cores=NCORES, overwrite_data=False)
