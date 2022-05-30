import os

import numpy as np


import analyse_functions_neil as af 

from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

base_dir = os.path.dirname(os.path.realpath(__file__))

cb = IscaCodeBase.from_directory(GFDL_BASE)


# expname 
#expname_base = 'GREY_T42L30fr_Drnl_om_1_2_a_1_1_ps_1_1_TJ001' # next 9/10
expname = lambda om_fact:'VENUS_T21L52_olddiff_highsponge_ures'
timestep  = 120.#60.#120

# planet params 
diurnal_period = -116.75
solar_constant = 150.*4
del_sol        = 1.#0.9#1.4
omega          = 2*np.pi / 86400. / 243#(-1*diurnal_period)
radius         = 6.051e6 #/ 20.
ps             = 92.e5 #*10.#
grav           = 8.87
R_d            = 188. 
cp             = 900.

# BL params 
gust        = 0.0 
rough_mom   = 5.e-3
rough_moist = 1.e-5
rough_sens  = 1.e-5
background_m = 0.01#5
background_t = 0.01#5
tj_bl_pres     = 8.5e4 * ps / 1.e5 
tj_strato_pres = 1.e4 * ps / 1.e5 
tj_coeff       = 0.001#2
sigma_b        = 0.7 
kf_days        = 1.

# sponge params
sponget = 0.5 
spongep = ps * 1.3e-5 #/ 20.

# ml params 
mldepth = 2.5 
mlalbedo = 0.0

# conv params 
gamma_conv = .77
dt_conv = 7200. 
do_relax = False

## RT params
tau_zero                = 1275#5.#4.54#4.47 
pres_zero               = ps##150000
do_rc_optical_depth     = False 
tau_rc                  = -1. 
pres_rc                 = -1.
linear_tau              = 0.0
nonlin_tau_exponent     = 1.#4./3.
lin_tau_exponent        = 1.

sw_k1                   = 0.0#141.22/4.#35.167 
sw_k2                   = 0.0#0.0643#0.0613 
do_log_sw = True
logf = 1 - 0.7
logtau_B  = 0.11204
logtau_sw_s = 0.90057
logtau_0 = 903.88
logtau_siglev = 0.005
logtau_pres_zero = ps
sw_tau_exponent = lin_tau_exponent

# these don't change 
bond_albedo             = 0.0#27 
sw_k1_frac              = 1.#4./7. #1.
do_sw_split             = False#True # False
del_sw                  = 0.0 
do_seasonal             = False 
solday                  = -10 
equinox_day             = 0 
dt_rad_avg              = -1 
use_time_average_coszen = False 
diabatic_acce           = 1.
tidally_locked          = False
noon_longitude          = 180.
simple_diurnal          = True
perp_eq                 = False


# core params 
fourier   = 21
spherical = fourier+1
lats      = 32
lons      = lats * 2
cores     = 16
levels    = 52
T_ini     = 500. 

# levels params 
# pk = np.array([40.00000,     100.00000,     200.00000,     
#               350.00000,     550.00000,     800.00000,     
#              1085.00000,    1390.00000,    1720.00000,     
#              2080.00000,    2470.00000,    2895.00000,     
#              3365.00000,    3890.00000,    4475.00000,     
#              5120.00000,    5830.00000,    6608.00000,     
#              7461.00000,    8395.00000,    9424.46289,     
#             10574.46900,   11864.80330,   13312.58850,     
#             14937.03770,   16759.70760,   18804.78670,     
#             21099.41250,   23674.03720,   26562.82650,     
#             29804.11680,   32627.31601,   34245.89759,     
#             34722.29104,   34155.20062,   32636.50533,     
#             30241.08406,   27101.45052,   23362.20912,     
#             19317.04955,   15446.17194,   12197.45091,     
#              9496.39912,    7205.66920,    5144.64339,     
#              3240.79521,    1518.62245,       0.00000,     
#                 0.00000]) * ps / 1.e5 
# bk = np.array([0.00000,       0.00000,       0.00000,      
#                0.00000,       0.00000,       0.00000,      
#                0.00000,       0.00000,       0.00000,      
#                0.00000,       0.00000,       0.00000,      
#                0.00000,       0.00000,       0.00000,      
#                0.00000,       0.00000,       0.00000,      
#                0.00000,       0.00000,       0.00000,      
#                0.00000,       0.00000,       0.00000,      
#                0.00000,       0.00000,       0.00000,      
#                0.00000,       0.00000,       0.00000,      
#                0.00000,       0.00813,       0.03224,      
#                0.07128,       0.12445,       0.19063,      
#                0.26929,       0.35799,       0.45438,      
#                0.55263,       0.64304,       0.71703,      
#                0.77754,       0.82827,       0.87352,      
#                0.91502,       0.95235,       0.98511,      

#                1.00000  ])

#ztwid = np.flip(np.arange(levels+1))/levels
#bk = np.exp(-5.*(0.05*ztwid + 0.95*ztwid**3.))
#pk = np.zeros_like(bk)
#new_pk = pk.tolist() #(pk + ps * bk).tolist()
#new_bk = bk.tolist() #(np.zeros_like(bk)).tolist()


cb.compile()  # compile the source code to working directory $GFDL_WORK/codebase

# create an Experiment object to handle the configuration of model parameters
# and output diagnostics
exp = Experiment(expname(1.), codebase=cb)
exp.inputfiles = [os.path.join(base_dir, 'create_ic/venus_col_log_check_newFlux_52H145_icwithu_T21.nc')]
#exp.inputfiles = [os.path.join(base_dir, 'create_ic/titan_col_log_check_newFlux_ic.nc')]

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 360, 'days', time_units='days')
#diag.add_file('atmos_monthly', 8, 'hours', time_units='days')

#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'height', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
#diag.add_field('dynamics', 'omega', time_avg=True)
diag.add_field('titan_gray_rad', 'olr', time_avg=True)
diag.add_field('titan_gray_rad', 'swdn_toa', time_avg=True)
# diag.add_field('titan_gray_rad', 'flux_lw', time_avg=True)
# diag.add_field('titan_gray_rad', 'flux_sw', time_avg=True)
# diag.add_field('titan_gray_rad', 'net_lw_surf', time_avg=True)
# diag.add_field('titan_gray_rad', 'swdn_sfc', time_avg=True)
# diag.add_field('mixed_layer', 'flux_t', time_avg=True)
# diag.add_field('atmosphere','dt_tg_convection',time_avg=True)
# diag.add_field('atmosphere','dt_tg_diffusion',time_avg=True)
# diag.add_field('rayleigh_bottom_drag','udt_rd',time_avg=True)
# diag.add_field('rayleigh_bottom_drag','vdt_rd',time_avg=True)

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
     'dt_atmos':timestep,
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
        'roughness_mom':rough_mom,
        'roughness_heat':rough_sens,
        'roughness_moist':rough_moist,                
        'two_stream_gray': False,     
        'titan_gray':True, 
        'convection_scheme': 'DRYADJ'
    },



    'vert_turb_driver_nml': {
        'do_mellor_yamada': False,     
        'do_diffusivity': True,        
        'do_edt'        :False, 
        'do_simple'    :  True,         
        'constant_gust' : gust,          
        'use_tau'       : False
    },

    'constants_nml':{
        # set planet constants 
        'omega'   : omega,
        'grav'    : grav, 
        'radius'  : radius, 
        'rdgas'   : R_d, 
        'PSTD'    : ps * 10., 
        'PSTD_MKS': ps, 
        'kappa'   : R_d / cp,
        'es0'     : 1.e-6, # small number just in case moisture tries to creep into model 
    },
    
    'diffusivity_nml': {
        'do_entrain':False,
        'do_simple' :False,
        'tj_bl_pres':tj_bl_pres,
        'tj_strato_pres':tj_strato_pres,
        'tj_coeff':tj_coeff,
        'background_m':background_m, 
        'background_t':background_t
    },

    'surface_flux_nml': {
        'use_virtual_temp': False,
        'do_simple' : True,
        'old_dtaudv': True, 
        'gust_const': gust,
        'tj_coeff':tj_coeff
    },

    'rayleigh_bottom_drag_nml':{
       'sigma_b':sigma_b,
       'kf_days':kf_days
    },
    

    'atmosphere_nml': {
        'idealized_moist_model': True
    },


    'mixed_layer_nml': {
        'tconst'                : T_ini+1.,
        'prescribe_initial_dist':True,
        'delta_T':0.0,
        'evaporation'           :False,   
        'depth'                 : mldepth,                  
        'albedo_value'          : mlalbedo,                    
    },

    'dry_adj_nml': {
        'alpha':gamma_conv,
        'tau_relax':dt_conv, 
        'do_relax':do_relax
    },

    
    'lscale_cond_nml': {
        'do_simple':True,
        'do_evap':True
    },
    
    'sat_vapor_pres_nml': {
        'do_simple':True,
        'tcmin_simple':-273.,
        'tcmax_simple':1000.,
    },
    
    'damping_driver_nml': {
        'do_rayleigh': True,
        'trayfric': -1*sponget,              # neg. value: time in *days*
        'sponge_pbottom':  spongep,           #(units are Pa)
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
        'lin_tau_exponent':lin_tau_exponent,
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
        'logf':logf,
        'sw_tau_exponent':sw_tau_exponent,
        'logtau_B': logtau_B, 
        'logtau_sw_s':logtau_sw_s, 
        'logtau_0': logtau_0, 
        'logtau_siglev':logtau_siglev,
        'logtau_pres_zero': logtau_pres_zero,
        'tidally_locked':tidally_locked, 
        'simple_diurnal':simple_diurnal,
        'diurnal_period':diurnal_period,
        'noon_longitude_tl':noon_longitude,
        'perp_eq':perp_eq
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
        'damping_coeff':1/86400.,   
        'do_water_correction':False, 
        'reference_sea_level_press':ps,
        'num_levels':levels,               #How many model pressure levels to use
        'valid_range_t':[1.,1000.],
        'initial_sphum':[0.0],#2.e-6],
        'vert_coord_option':'hybrid', #Use the vertical levels from Hammond and Lewis 2021 
        'robert_coeff':0.03, 
        'lon_max'                  : lons,
        'lat_max'                  : lats,
        'num_fourier'              : fourier,
        'num_spherical'            : spherical,
        'surf_res':0.15,
        'scale_heights' : 14.5,
        'exponent':3.0,
    },

        'spectral_init_cond_nml': {
        'initial_temperature': T_ini, # initial atmospheric temperature 
    },

    
    'ic_from_external_file_nml':{ 
        'file_name':'INPUT/venus_col_log_check_newFlux_52H145_icwithu_T21.nc', 
        'u_name': 'ucomp', 
        'v_name': 'vcomp', 
        't_name': 'temp', 
        'ps_name': 'ps'
    },
    
    #'vert_coordinate_nml': {
     #    'bk': new_bk,
     #    'pk': new_pk,
     #}
})

#Lets do a run!
if __name__=="__main__":

    run_exp = exp.derive(expname(1.))
    run_exp.namelist['spectral_dynamics_nml']['initial_state_option'] = 'input'
    ds = af.open_experiment('venus_col_log_check_newFlux_52H145',50,50,'atmos_monthly.nc',
                            variables=['temp'],cal=False)
    T_ini_sfc = np.squeeze(ds.temp.values[-1,-1])
    run_exp.namelist['mixed_layer_nml']['tconst'] = T_ini_sfc+1.
    run_exp.run(1, num_cores=cores, use_restart=False, overwrite_data=False)
    for i in range(2,201):
        exp.run(i, num_cores=cores, overwrite_data=False)
            


