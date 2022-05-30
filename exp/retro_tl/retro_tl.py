import os

import numpy as np

from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE, FailedRunError

NCORES = 16
base_dir = os.path.dirname(os.path.realpath(__file__))
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
exp = Experiment('retro_tl_150', codebase=cb)
#exp = Experiment('byrne_sc_test', codebase=cb)



#Empty the run directory ready to run
exp.clear_rundir()


#Define values for the 'core' namelist
exp.namelist = namelist = Namelist({
    'main_nml':{
     'days': 25,
     'dt_atmos': 600.,
     'calendar' : 'no_calendar',
     'current_date' : [1,1,1,0,0,0]
    },

    'idealized_moist_phys_nml': {
        'do_damping': True,
        'turb':True,
        'mixed_layer_bc':True,
        'do_virtual' :False,
        'do_simple': True,
        'roughness_mom':3.21e-05,
        'roughness_heat':3.21e-05,
        'roughness_moist':3.21e-05,                
        'two_stream_gray': True,     #Use grey radiation
	    'do_rrtm_radiation': False, 
        'convection_scheme': 'dry', #Use the simple Betts Miller convection scheme from Frierson
    },

    'vert_turb_driver_nml': {
        'do_mellor_yamada': False,     # default: True
        'do_diffusivity': True,        # default: False
        'do_simple': True,             # default: False
        'constant_gust': 0.0,          # default: 1.0
        'use_tau': False
    },
    
    'diffusivity_nml': {
        'do_entrain':False,
        'do_simple': True,
    },


    'two_stream_gray_rad_nml': {
        'rad_scheme': 'byrne',    
        'bog_a' :1.0, 
        'bog_b': 0.0,
        'do_seasonal': True,                
        'atm_abs': 0.0,  
        'solar_constant':1360.0,                         
    },

    'surface_flux_nml': {
        'use_virtual_temp': False,
        'do_simple': False,
        'old_dtaudv': True    
    },

    'atmosphere_nml': {
        'idealized_moist_model': True
    },

    'mixed_layer_nml': {
        'tconst' : 280.,
        'prescribe_initial_dist':True,
        'evaporation':False,   
        'depth': 5.0,                          #Depth of mixed layer used
        'albedo_value': 0.3,                  #Albedo value used             
        'do_qflux': False,
    },

    'dry_convection_nml': {
        'tau':7200.,
        'Gamma':1.0, 

    },

    
    
    'damping_driver_nml': {
        'do_rayleigh': True,
        'trayfric': -0.25,              # neg. value: time in *days*
        'sponge_pbottom':  5000.,           #Bottom of the model's sponge down to 50hPa (units are Pa)
        'do_conserve_energy': True,             
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

    'sat_vapor_pres_nml': { 
        'do_simple': True, 
        'tcmin_simple': -200, 
        'tcmax_simple': 300
    },

    'spectral_dynamics_nml': {
        'damping_order': 4,  
        'reference_sea_level_press':1.0e5,
        'num_levels':25,              
        'valid_range_t':[40.,800.],
        'initial_sphum':[0.0],
        'vert_coord_option':'uneven_sigma', 
        'surf_res':0.2,
        'scale_heights' : 11.0,
        'exponent':7.0,
        'robert_coeff':0.03, 
        'do_water_correction': False
    },

    'astronomy_nml': {
        'ecc': 0.0,
        'obliq': 0.0,     # planetary obliquity
        'per': 0.0,        # begin at autumn equinox

    },
    
    'constants_nml': { 
        'orbital_period': 10. * 86400., # 10 days 
        'omega': 2.*3.1415926/(10.*86400.), 
        'rdgas': 297.,
        'kappa': 297. / 1039.0,
        'wtmair': 28.013
        #'CP_AIR': 900.,

    },
})

#Lets do a run!
if __name__=="__main__":
    resolution = 'T42', 25
    exp.set_resolution(*resolution)
    #Tell model how to write diagnostics
    diag = DiagTable()
    diag.add_file('atmos_monthly', 5, 'days', time_units='days')
    #Tell model which diagnostics to write
    diag.add_field('dynamics', 'ps', time_avg=True)
    diag.add_field('atmosphere', 'rh', time_avg=True)
    diag.add_field('dynamics', 'bk')
    diag.add_field('dynamics', 'pk')
    diag.add_field('mixed_layer', 't_surf', time_avg=True)
    diag.add_field('mixed_layer', 'flux_lhe', time_avg=True)
    diag.add_field('mixed_layer', 'flux_t', time_avg=True)
    diag.add_field('dynamics', 'sphum', time_avg=True)
    diag.add_field('dynamics', 'ucomp', time_avg=True)
    diag.add_field('dynamics', 'vcomp', time_avg=True)
    diag.add_field('dynamics', 'omega', time_avg=True)
    diag.add_field('dynamics', 'temp', time_avg=True)
    diag.add_field('dynamics', 'vor', time_avg=True)
    diag.add_field('dynamics', 'div', time_avg=True)
    diag.add_field('dynamics', 'height', time_avg=True)
    diag.add_field('dynamics', 'height_half', time_avg=True)
    diag.add_field('dynamics', 'pres_full', time_avg=True)
    diag.add_field('dynamics', 'pres_half', time_avg=True)
    diag.add_field('two_stream', 'swdn_toa', time_avg=True)
    exp.diag_table = diag
    
    exp.run(1, use_restart=False, num_cores=NCORES, overwrite_data=False)      
    for i in range(2,41):   
        exp.run(i, num_cores=NCORES, overwrite_data=False)