import os

import numpy as np

from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES =32

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
exp = Experiment('tls_test_experiment', codebase=cb)

#exp.inputfiles = [os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc')]

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_yearly', 360, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)



exp.diag_table = diag


#Empty the run directory ready to run
exp.clear_rundir()

#Define values for the 'core' namelist
exp.namelist = namelist = Namelist({
    'main_nml': {
        'days'   : 360,
        'hours'  : 0,
        'minutes': 0,
        'seconds': 0,
        'dt_atmos':900,
        'current_date' : [1,1,1,0,0,0],
        'calendar' : 'thirty_day'
    },

    'idealized_moist_phys_nml': {
        'two_stream_gray': False,
        'do_rrtm_radiation': True,    #Use RRTM radiation, not grey
        'convection_scheme': 'SIMPLE_BETTS_MILLER',     #Use the simple Betts Miller convection scheme
        'do_damping': True,
        'turb':True,
        'mixed_layer_bc':True,
        'do_virtual' :False,
        'do_simple': True,
        'roughness_mom':3.21e-05,
        'roughness_heat':3.21e-05,
        'roughness_moist':3.21e-05,                
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
        'do_simple': False,
    },

    'surface_flux_nml': {
        'use_virtual_temp': False,
        'do_simple': True,
        'old_dtaudv': True    
    },

    'atmosphere_nml': {
        'idealized_moist_model': True
    },

    #Use a large mixed-layer depth, and the Albedo of the CTRL case in Jucker & Gerber, 2017
    'mixed_layer_nml': {
        'depth': 30,
        'albedo_value': 0.25,
        'tconst' : 285.,
        'prescribe_initial_dist':True,
        'evaporation':True,
        'do_qflux': False,     
    },

    'qe_moist_convection_nml': {
        'rhbm':0.7,
        'Tmin':160.,
        'Tmax':350.   
    },
    
    'lscale_cond_nml': {
        'do_simple':True,
        'do_evap':True
    },
    
    'sat_vapor_pres_nml': {
        'do_simple':True
    },
    
    'damping_driver_nml': {
        'do_rayleigh': True,
        'trayfric': -0.5,              # neg. value: time in *days*
        'sponge_pbottom':  50.,
        'do_conserve_energy': True,         
    },

#     'qflux_nml': {
#         'qflux_amp': 30.0
#     },

    'rrtm_radiation_nml': {
        'solr_cnst': 1360.,  #s set solar constant to 1360, rather than default of 1368.22
        'dt_rad': 7200, #Use long RRTM timestep
        'do_read_ozone':False,
        'frierson_solar_rad':True, 
        'del_sol':1.2,
        'include_secondary_gases':True,
        'co2ppmv':355.,
        'n2o_val':320.e-9, 
        'ch4_val':1700.e-9, 
        'o3_val':30.e-9, # NEED TO ADD THIS 
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
        'num_levels':40,
        'valid_range_t':[100.,800.],
        'initial_sphum':[2.e-6],
        'vert_coord_option':'uneven_sigma',
        'surf_res':0.5,
        'scale_heights' : 11.0,
        'exponent':7.0,
        'robert_coeff':0.03
    }
    
    
})

# Lets do a run!
if __name__=="__main__":
    
    run_exp = exp.derive('tls_ctrl')
    
    run_exp.run(1, use_restart=False, num_cores=NCORES, overwrite_data=False)
    for i in range(2,31):
        if i == 21:
            diag = DiagTable()
            diag.add_file('atmos_6hrly', 6, 'hours', time_units='days')

            #Write out diagnostics need for vertical interpolation post-processing
            diag.add_field('dynamics', 'ps', time_avg=True)
            diag.add_field('dynamics', 'bk')
            diag.add_field('dynamics', 'pk')
            diag.add_field('dynamics', 'zsurf')

            #Tell model which diagnostics to write
            diag.add_field('mixed_layer', 't_surf', time_avg=True)
            diag.add_field('dynamics', 'temp', time_avg=True)
            diag.add_field('atmosphere',   'temp_2m', time_avg=True)
            diag.add_field('atmosphere',   'sphum_2m', time_avg=True)
            diag.add_field('atmosphere',   'u_10m', time_avg=True)
            diag.add_field('atmosphere',   'v_10m', time_avg=True) 

            run_exp.diag_table = diag 
        run_exp.run(i, num_cores=NCORES, overwrite_data=False)
        
    run_exp2 = exp.derive('tls_4xCO2')
    run_exp2.namelist['rrtm_radiation_nml']['co2ppmv'] = 355*4.
    
    run_exp2.run(1, use_restart=False, num_cores=NCORES, overwrite_data=False)
    for i in range(2,31):
        if i == 21:
            diag = DiagTable()
            diag.add_file('atmos_6hrly', 6, 'hours', time_units='days')

            #Write out diagnostics need for vertical interpolation post-processing
            diag.add_field('dynamics', 'ps', time_avg=True)
            diag.add_field('dynamics', 'bk')
            diag.add_field('dynamics', 'pk')
            diag.add_field('dynamics', 'zsurf')

            #Tell model which diagnostics to write
            diag.add_field('mixed_layer', 't_surf', time_avg=True)
            diag.add_field('dynamics', 'temp', time_avg=True)
            diag.add_field('atmosphere',   'temp_2m', time_avg=True)
            diag.add_field('atmosphere',   'sphum_2m', time_avg=True)
            diag.add_field('atmosphere',   'u_10m', time_avg=True)
            diag.add_field('atmosphere',   'v_10m', time_avg=True) 

            run_exp2.diag_table = diag 
        run_exp2.run(i, num_cores=NCORES, overwrite_data=False)
