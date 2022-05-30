import os

import numpy as np

from isca import GreyCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES = 32
base_dir = os.path.dirname(os.path.realpath(__file__))
# a CodeBase can be a directory on the computer,
# useful for iterative development
cb = GreyCodeBase.from_directory(GFDL_BASE)

# or it can point to a specific git repo and commit id.
# This method should ensure future, independent, reproducibility of results.
# cb = DryCodeBase.from_repo(repo='https://github.com/isca/isca', commit='isca1.1')

# compilation depends on computer specific settings.  The $GFDL_ENV
# environment variable is used to determine which `$GFDL_BASE/src/extra/env` file
# is used to load the correct compilers.  The env file is always loaded from
# $GFDL_BASE and not the checked out git repo.

gust = 0.0

## RT params 
tau_zero                = 0.5#15
pres_zero               = 1.e5
do_rc_optical_depth     = False
tau_rc                  = (0.11/1.15)*tau_zero
pres_rc                 = (30403./1.e5)*pres_zero
linear_tau              = 1.0
nonlin_tau_exponent     = 1.
sw_k1                   = 0.0 #51.8 #35.168 #0.0 #
sw_k2                   = 0.0 #0.071 #0.0599 # #
sw_k3                   = 0.0 #8.0745  / 1.305e2 * sw_k1 
sw_gamma                = 3.395e-1

# these don't change 
bond_albedo             = 0.0#27 
sw_k1_frac              = 4./7. 
do_sw_split             = True 
solar_constant          = 1000.
del_sol                 = 1.2
del_sw                  = 0.0 
do_seasonal             = False
solday                  = -10 
equinox_day             = 0 
dt_rad_avg              = -1 
use_time_average_coszen = False 
diabatic_acce           = 1.
tidally_locked          = True 
noon_longitude          = 180.

# levels params 


pk = np.array([40.00000,     100.00000,     200.00000,     
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
                0.00000]) 
bk = np.array([0.00000,       0.00000,       0.00000,      
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
               1.00000  ])
nlevs = 48

# bk = np.array([0.000000, 0.0117665, 0.0196679, 0.0315244, 0.0485411, 0.0719344, 0.1027829, 0.1418581, 0.1894648, 0.2453219, 0.3085103, 0.3775033, 0.4502789, 0.5244989, 0.5977253, 0.6676441, 0.7322627, 0.7900587, 0.8400683, 0.8819111, 0.9157609, 0.9422770, 0.9625127, 0.9778177, 0.9897489, 1.0000000]) 
# pk = np.zeros_like(bk)
# nlevs = 25

new_pk = pk.tolist() #(pk + ps * bk).tolist()
new_bk = bk.tolist() #(np.zeros_like(bk)).tolist()


cb.compile()  # compile the source code to working directory $GFDL_WORK/codebase

# create an Experiment object to handle the configuration of model parameters
# and output diagnostics
exp = Experiment('tl_test_experiment_tj_0000001', codebase=cb)

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('mixed_layer', 't_surf', time_avg=True)
#diag.add_field('mixed_layer', 'flux_t', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'height', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'omega', time_avg=True)
#diag.add_field('dynamics', 'dsefluxdiv', time_avg=True)
#diag.add_field('dynamics', 'dsefluxdivint', time_avg=True)
#diag.add_field('dynamics', 'rt_dlnpsdt_int', time_avg=True)
#diag.add_field('dynamics', 'dt_tg_core', time_avg=True)
#diag.add_field('dynamics', 'dt_tg_damp', time_avg=True)
#diag.add_field('atmosphere', 'dt_tg_diffusion', time_avg=True)
#diag.add_field('atmosphere', 'dt_tg_phys', time_avg=True)
diag.add_field('titan_gray_rad', 'olr', time_avg=True)
diag.add_field('titan_gray_rad', 'swdn_toa', time_avg=True)
#diag.add_field('titan_gray_rad', 'swdn_sfc', time_avg=True)
#diag.add_field('titan_gray_rad', 'net_lw_surf', time_avg=True)
#diag.add_field('titan_gray_rad', 'tdt_rad', time_avg=True)


exp.diag_table = diag

#Empty the run directory ready to run
exp.clear_rundir()

#Define values for the 'core' namelist
exp.namelist = namelist = Namelist({
    'main_nml':{
     'days'   : 360,
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,
     'dt_atmos':600,
     'current_date' : [1,1,1,0,0,0],
     'calendar' : 'thirty_day'
    },

    'idealized_moist_phys_nml': {
        'do_damping': True,
        'turb':True,
        'mixed_layer_bc':True,
        'do_tj_bl':True,
        'terr_ray_surf':True, 
        'do_virtual' :False,
        'do_simple': True,
        'roughness_mom':3.21e-05,
        'roughness_heat':3.21e-05,
        'roughness_moist':3.21e-05,                
        'two_stream_gray': False,     #Use grey radiation
        'titan_gray':True, 
        'convection_scheme': 'NONE', #Use the simple Betts Miller convection scheme from Frierson
    },

    'rayleigh_bottom_drag_nml':{ 
        'sigma_b':0.7, 
        'kf_days':1.0, 

    },

    'vert_turb_driver_nml': {
        'do_mellor_yamada':False,     # default: True
        'do_diffusivity': False,        # default: False
        'do_edt':False, 
        'do_simple': True,             # default: False
        'constant_gust': gust,          # default: 1.0
        'use_tau': True
    },

    'constants_nml':{ 
        'es0':0.0, 
        'omega':7.29e-5 / 10.
    }, 
    
    'diffusivity_nml': {
        'do_entrain':False,
        'do_simple': True,
        'fixed_depth':False, 
        'rich_crit_pbl':1.0,#0.1
        'tj_bl_pres':8.5e4, 
        'tj_strato_pres':1.e4,
        'tj_coeff': 0.0000001
        #'depth_0':8000.0
    },

    'surface_flux_nml': {
        'use_virtual_temp': False,
        'do_simple': True,
        'old_dtaudv': True, 
        'gust_const': gust, 
        'tj_coeff': 0.0000001
    },

    'atmosphere_nml': {
        'idealized_moist_model': True
    },

    #Use a large mixed-layer depth, and the Albedo of the CTRL case in Jucker & Gerber, 2017
    'mixed_layer_nml': {
        'tconst' : 285.,
        'prescribe_initial_dist':True,
        'evaporation':True,   
        'depth': 2.5,                          #Depth of mixed layer used
        'albedo_value': 0.0,                  #Albedo value used             
    },

    'qe_moist_convection_nml': {
        'rhbm':0.7,
        'Tmin':160.,
        'Tmax':350.   
    },

    'betts_miller_nml': {
       'rhbm': .7   , 
       'do_simp': False, 
       'do_shallower': True, 
    },
    
    'lscale_cond_nml': {
        'do_simple':True,
        'do_evap':True
    },

    'llcs_nml': {
        'alpha':1.,
	    'llcs_timescale':7200.,
        'do_simple':True
    },
    
    'sat_vapor_pres_nml': {
        'do_simple':True,                        
        'do_not_calculate':False, 
        'tcmin_simple' : -99999,
    	'tcmax_simple' : 99999
    },
    
    'damping_driver_nml': {
        'do_rayleigh': True,
        'trayfric': -0.5,              # neg. value: time in *days*
        'sponge_pbottom':  1000.,           #Bottom of the model's sponge down to 50hPa (units are Pa)
        'do_conserve_energy': True,             
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
        'tidally_locked':tidally_locked, 
        'noon_longitude_tl':noon_longitude,
        'simple_diurnal':True
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

    'spectral_dynamics_nml': {
        'damping_order': 4,             
        'water_correction_limit': 200.e2,
        'reference_sea_level_press':1.0e5,
        'num_levels':nlevs,               #How many model pressure levels to use
        'valid_range_t':[100.,800.],
        'initial_sphum':[0.0],#[2.e-6],
        'vert_coord_option':'input', #Use the vertical levels from Frierson 2006
        #'surf_res':0.5,
        #'scale_heights' : 11.0,
        #'exponent':7.0,
        'robert_coeff':0.03, 
        'lon_max'                  : 128,
        'lat_max'                  : 64,
        'num_fourier'              : 42,
        'num_spherical'            : 43,
        #'vert_advect_t':'FOURTH_CENTERED',
        #'vert_advect_uv':'FOURTH_CENTERED'
    },
    'vert_coordinate_nml': {
        'bk': new_bk,
        'pk': new_pk,
       }
})

#Lets do a run!
if __name__=="__main__":
    exp.run(1, use_restart=False, num_cores=NCORES, overwrite_data=True)
    for i in range(2,6):#121):
        exp.run(i, num_cores=NCORES, overwrite_data=True)
