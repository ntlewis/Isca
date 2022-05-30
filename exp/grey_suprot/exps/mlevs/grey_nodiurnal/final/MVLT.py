import numpy as np

from isca import GreyCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

cb = GreyCodeBase.from_directory(GFDL_BASE)

#expname
expname = 'mv_T42L20_om_1_1_a_1_20_ps_1_1'
timestep  = 60.#120#60#60.#120

# planet params 
omega          = 7.27e-5 
radius         = 6.371e6 /20.#/ 20.
ps             = 1.e5 
grav           = 9.81
R_d            = 287. 
cp             = 1005.

# BL params 
gust        = 0.0 
rough_mom   = 3.21e-5
rough_moist = 3.21e-5
rough_sens  = 3.21e-5

# sponge params
sponget = 0.5 
spongep = 5000. 

# ml params 
mldepth = 2.5 
mlalbedo = 0.0


# core params 
fourier   = 42 
spherical = fourier+1
lats      = 64
lons      = lats * 2
cores     = 32
levels    = 20#48
T_ini     = 280. 

cb.compile()  

exp = Experiment(expname, codebase=cb)

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 1, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'omega', time_avg=True)
diag.add_field('dynamics', 'height', time_avg=True)
diag.add_field('hs_forcing', 'tdt', time_avg=True)
diag.add_field('hs_forcing', 'udt_rdamp', time_avg=True)
diag.add_field('hs_forcing', 'vdt_rdamp', time_avg=True)

exp.diag_table = diag

#Empty the run directory ready to run
exp.clear_rundir()

# define namelist values as python dictionary
# wrapped as a namelist object.
namelist = Namelist({
    'main_nml': {
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
        'do_tj_bl':False,
        'terr_ray_surf':False, 
        'do_virtual' :False,
        'do_simple': True,
        'roughness_mom':rough_mom,
        'roughness_heat':rough_sens,
        'roughness_moist':rough_moist,             
        'two_stream_gray': False,
        'titan_gray':False, 
        'do_newtonian_cooling':True, 
        'convection_scheme': 'NONE', 
    },

    'vert_turb_driver_nml': {
        'do_mellor_yamada': False,     # default: True
        'do_diffusivity': False,        # default: False
        'do_simple': True,             # default: False
        'constant_gust': gust,          # default: 1.0
        'use_tau': False, 
        'do_background_only':False, 
        'diff_m_background_only':0.0,#0.01, 
        'diff_t_background_only':0.0,#0.01,
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
        'do_simple' : True,
    },

    'surface_flux_nml': {
        'use_virtual_temp': False,
        'do_simple' : True,
        'old_dtaudv': True, 
        'gust_const': gust
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

    'lscale_cond_nml': {
        'do_simple':True,
        'do_evap':True
    },
    
    'sat_vapor_pres_nml': {
        'do_simple':True,   
        'do_not_calculate':False, 
        'tcmin_simple' : -99999,
    	'tcmax_simple' : 99999
    },
    
    'damping_driver_nml': {
        'do_rayleigh': True,
        'trayfric': -1*sponget,              # neg. value: time in *days*
        'sponge_pbottom':  spongep,           #(units are Pa)
        'do_conserve_energy': True,             
    },

    # configure the relaxation profile
    'hs_forcing_nml':{
        'equilibrium_t_option':'Mitchell_Vallis', #forcing like Mitchell and Vallis (2010), J. Geophys. Res. 
        't_zero': 285.,    # mean surface temp
        'delh': 0.2,       # non-dim equator-pole temp gradient 
        'alpha': 0.176,    # sets R/cp used by lapse rate > this yields moist adiabat 
        'sigma_b': 0.7,    # boundary layer friction height (default p/ps = sigma = 0.7)

        # negative sign is a flag indicating that the units are days
        'ka':   -40.,      # Constant Newtonian cooling timescale (default 40 days)
        'ks':    -4.,      # Boundary layer dependent cooling timescale (default 4 days)
        'kf':   -1.,       # BL momentum frictional timescale (default 1 days)

        'do_conserve_energy':   True,  # convert dissipated momentum into heat (default True)
    },

    
    
    'spectral_dynamics_nml': {
        'damping_order': 4,             
        'do_water_correction':False, 
        'reference_sea_level_press':ps,
        'num_levels':levels,               #How many model pressure levels to use
        'valid_range_t':[100.,800.],
        'initial_sphum':[0.0],#2.e-6],
        'vert_coord_option':'even_sigma', #Use the vertical levels from Hammond and Lewis 2021 
        'robert_coeff':0.03, 
        'lon_max'                  : lons,
        'lat_max'                  : lats,
        'num_fourier'              : fourier,
        'num_spherical'            : spherical,
    },

        'spectral_init_cond_nml': {
        'initial_temperature': T_ini, # initial atmospheric temperature 
    },


    'diag_manager_nml': {
        'mix_snapshot_average_fields': False
    },

    'fms_nml': {
        'domains_stack_size': 600000                        # default: 0
    },

    'fms_io_nml': {
        'threading_write': 'single',                         # default: multi
        'fileset_write': 'single',                           # default: multi
    }
})

exp.namelist = namelist

#Lets do a run!
if __name__ == '__main__':
    exp.run(1, num_cores=cores, use_restart=False, overwrite_data=False)
    for i in range(2, 11):
        exp.run(i, num_cores=cores, overwrite_data=False)  