import os

import numpy as np

from isca import GreyCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE



cb = GreyCodeBase.from_directory(GFDL_BASE)


# expname 
#expname_base = 'GREY_T42L30fr_Drnl_om_1_2_a_1_1_ps_1_1_TJ001' # next 9/10
expname = lambda om_fact:'GREY_T42L40fr6_Drnl_om_1_'+str(int(om_fact))+'_a_1_1_ps_10_1_TJ001_n1_abs1000hPahalffrac'
timestep  = 360.#60.#120

# planet params 
diurnal_period = -10.
solar_constant = 1000. 
del_sol        = 1.#0.9#1.4
omega          = -2*np.pi / 86400. / diurnal_period
radius         = 6.371e6 #/ 20.
ps             = 1.e6 #*10.#
grav           = 9.81
R_d            = 287. 
cp             = 1005.

# BL params 
gust        = 0.0 
rough_mom   = 5.e-3
rough_moist = 1.e-5
rough_sens  = 1.e-5
tj_bl_pres     = 8.5e4 * ps / 1.e5 
tj_strato_pres = 1.e4 * ps / 1.e5 
tj_coeff       = 0.001
sigma_b        = 0.7 
kf_days        = 1.

# sponge params
sponget = 0.5 
spongep = ps * 0.01 #/ 20.

# ml params 
mldepth = 2.5 
mlalbedo = 0.0

# conv params 
gamma_conv = 1.
dt_conv = 7200. 
do_relax = False

## RT params 
tau_zero                = 1.#0.125 #0.5#15
pres_zero               = ps#(31622./1.e5)*ps #1.e5
do_rc_optical_depth     = False
tau_rc                  = -1#(0.11/1.15)*tau_zero
pres_rc                 = -1#(30403./1.e5)*pres_zero
linear_tau              = 0.0#.25 #1.0
nonlin_tau_exponent     = 1.0 #4./3.#1.#4./3.
sw_k1                   = 1/(0.1**(1.))#1/(0.2**(4./3.)) #51.8 #35.168 #0.0 #
sw_k2                   = 0.0 #0.071 #0.0599 # #
sw_k3                   = 0.0 #8.0745  / 1.305e2 * sw_k1 
sw_gamma                = 0.0#3.395e-1

# these don't change 
bond_albedo             = 0.0#27 
sw_k1_frac              = 0.5#3./10. 
do_sw_split             = True 
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
fourier   = 42
spherical = fourier+1
lats      = 64
lons      = lats * 2
cores     = 32
levels    = 40
T_ini     = 280. 

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

ztwid = np.flip(np.arange(levels+1))/levels
bk = np.exp(-6.*(0.05*ztwid + 0.95*ztwid**3.))
pk = np.zeros_like(bk)
new_pk = pk.tolist() #(pk + ps * bk).tolist()
new_bk = bk.tolist() #(np.zeros_like(bk)).tolist()


cb.compile()  # compile the source code to working directory $GFDL_WORK/codebase

# create an Experiment object to handle the configuration of model parameters
# and output diagnostics
exp = Experiment(expname(1.), codebase=cb)

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
diag.add_field('dynamics', 'omega', time_avg=True)
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
        'tj_coeff':tj_coeff
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
        'tcmax_simple':800.,
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
        'do_water_correction':False, 
        'reference_sea_level_press':ps,
        'num_levels':levels,               #How many model pressure levels to use
        'valid_range_t':[100.,800.],
        'initial_sphum':[0.0],#2.e-6],
        'vert_coord_option':'input', #Use the vertical levels from Hammond and Lewis 2021 
        'robert_coeff':0.03, 
        'lon_max'                  : lons,
        'lat_max'                  : lats,
        'num_fourier'              : fourier,
        'num_spherical'            : spherical,
    },

        'spectral_init_cond_nml': {
        'initial_temperature': T_ini, # initial atmospheric temperature 
    },

    'vert_coordinate_nml': {
         'bk': new_bk,
         'pk': new_pk,
     }
})

#Lets do a run!
if __name__=="__main__":
    om_facts = np.array([10.,100.])#2.,4.,6.,8.,10.,20.,40.,60.,80.,100.,200.])
    for om_fact in om_facts:
        expname_use = expname(om_fact)
        run_exp = exp.derive(expname_use)#+'_12hrly')
        run_exp.namelist['constants_nml']['omega'] = 2*np.pi / 86400. / om_fact
        run_exp.namelist['titan_gray_rad_nml']['diurnal_period'] = -1 * om_fact
        run_exp.run(1, use_restart=False, num_cores=cores, overwrite_data=False)#, restart_file='/network/group/aopp/planetary/PLR018_LEWIS_VENUS2/isca_data/'+expname_use+'/restarts/res0010.tar.gz')
        for i in range(2,151):#121):
            run_exp.run(i, num_cores=cores, overwrite_data=False)