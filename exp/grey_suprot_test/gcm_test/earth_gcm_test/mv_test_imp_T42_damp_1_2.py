import numpy as np

from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES = 32
RESOLUTION = 'T42', 20  # T42 horizontal resolution, 25 levels in pressure

# a CodeBase can be a directory on the computer,
# useful for iterative development
cb = IscaCodeBase.from_directory(GFDL_BASE)

# or it can point to a specific git repo and commit id.
# This method should ensure future, independent, reproducibility of results.
# cb = DryCodeBase.from_repo(repo='https://github.com/isca/isca', commit='isca1.1')

# compilation depends on computer specific settings.  The $GFDL_ENV
# environment variable is used to determine which `$GFDL_BASE/src/extra/env` file
# is used to load the correct compilers.  The env file is always loaded from
# $GFDL_BASE and not the checked out git repo.

cb.compile()  # compile the source code to working directory $GFDL_WORK/codebase

# create an Experiment object to handle the configuration of model parameters
# and output diagnostics

exp_name = 'mv_T42_L20_dt120_even_rad_1_8_imp_damp_1_2'
exp = Experiment(exp_name, codebase=cb)

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'vor', time_avg=True)
diag.add_field('dynamics', 'div', time_avg=True)
diag.add_field('atmosphere','dt_ug_diffusion', time_avg=True)
diag.add_field('atmosphere','dt_vg_diffusion', time_avg=True)
diag.add_field('atmosphere','dt_tg_diffusion', time_avg=True)

exp.diag_table = diag

# define namelist values as python dictionary
# wrapped as a namelist object.
namelist = Namelist({
    'main_nml': {
        'dt_atmos': 120,
        'days': 30,
        'calendar': 'thirty_day',
        'current_date': [2000,1,1,0,0,0]
    },

    'atmosphere_nml': {
        'idealized_moist_model': True  # False for Newtonian Cooling.  True for Isca/Frierson
    },

    'spectral_dynamics_nml': {
        'damping_order'           : 4,                      # default: 2 
        'damping_coeff': 1 / (1./5. * 86400.), #'eddy_sponge_coeff' : 1/(1./10. * 86400.),
        'water_correction_limit'  : 200.e2,                 # default: 0
        'reference_sea_level_press': 1.0e5,                  # default: 101325
        'valid_range_t'           : [100., 800.],           # default: (100, 500)
        'initial_sphum'           : 0.0,                  # default: 0
        'vert_coord_option'       : 'even_sigma',         # default: 'even_sigma'
    },

    'idealized_moist_phys_nml': {
        'do_damping': False,
        'turb':True,
        'mixed_layer_bc':True,
        'do_virtual' :False,
        'do_simple': True,
        'roughness_mom':5.e-3,
        'roughness_heat':5.e-3,
        'roughness_moist':1.e-5,                
        'two_stream_gray': False,     #Use grey radiation
        'do_newtonian_cooling':True, 
        'convection_scheme': 'NONE', #Use the simple Betts Miller convection scheme from Frierson
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

    'vert_turb_driver_nml': {
        'do_mellor_yamada': False,     # default: True
        'do_diffusivity': False,        # default: False
        'do_simple': True,             # default: False
        'constant_gust': 0.0,          # default: 1.0
        'use_tau': False, 
        'do_background_only':False, 
        'diff_m_background_only':0.0,#0.01, 
        'diff_t_background_only':0.0,#0.01,
    },



        
    'damping_driver_nml': {
        'do_rayleigh': False,
        'trayfric': -(1./2.),              # neg. value: time in *days*
        'sponge_pbottom':  50.,           #Bottom of the model's sponge down to 50Pa (units are Pa)
        'do_conserve_energy': False,            
        'rayleigh_eddy_only':False, 
    },

    'lscale_cond_nml': {
        'do_simple':True,
        'do_evap': False 
    },
    
    'sat_vapor_pres_nml': {
        'do_simple':True,
        'tcmin_simple' : -273, 
        'tcmax_simple' : 2000
    },

    'mixed_layer_nml': {
        'tconst' : 285.,
        'prescribe_initial_dist':True,
        'evaporation':False,   
        'depth': 1.,                          #Depth of mixed layer used
        'albedo_value': 0.31,                  #Albedo value used             
    },



    'constants_nml':{
        'radius' : 6400.e3 / 8.
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
exp.set_resolution(*RESOLUTION)

#Lets do a run!
if __name__ == '__main__':
    exp.run(1, num_cores=NCORES, use_restart=False, overwrite_data=True)
    for i in range(2, 37):
        exp.run(i, num_cores=NCORES, overwrite_data=True)  # use the restart i-1 by default