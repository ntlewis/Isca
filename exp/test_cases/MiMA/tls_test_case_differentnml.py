import os

import numpy as np

from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES =8

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
        'do_virtual' :True,
        'do_simple': False,
        'roughness_mom':5.e-3,#3.21e-05,
        'roughness_heat':1.e-5,#3.21e-05,
        'roughness_moist':1.e-5,#3.21e-05,                
    },

    'vert_turb_driver_nml': {
        'do_mellor_yamada': False,     # default: True
        'do_diffusivity': True,        # default: False
        'do_simple': False,             # default: False
        'constant_gust': 1.0,          # default: 1.0
        'use_tau': False
    },
    
    'diffusivity_nml': {
        'do_entrain':False,
        'do_simple': False,
    },

    'surface_flux_nml': {
        'use_virtual_temp': True,
        'do_simple': False,
        'old_dtaudv': False    
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
        'tau_bm':7200.,
        'Tmin':120.,
        'Tmax':360.   
    },
    
    'lscale_cond_nml': {
        #'do_simple':True,
        'do_evap':False
    },
    
    'sat_vapor_pres_nml': {
        'do_simple':True
    },
    
    'damping_driver_nml': {
        'do_rayleigh': True,
        'trayfric': -0.5,              # neg. value: time in *days*
        'sponge_pbottom':  500.,
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
        'num_levels':48,
        'valid_range_t':[100.,800.],
        'initial_sphum':[2.e-6],
        'vert_coord_option':'input',
        #'surf_res':0.5,
        #'scale_heights' : 11.0,
        #'exponent':7.0,
        'robert_coeff':0.03
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
    
    
})

# Lets do a run!
if __name__=="__main__":
    
    run_exp = exp.derive('tls_diffnml_ctrl')
    
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
        
#     run_exp2 = exp.derive('tls_diffnml_4xCO2')
#     run_exp2.namelist['rrtm_radiation_nml']['co2ppmv'] = 355*4.
    
#     run_exp2.run(1, use_restart=False, num_cores=NCORES, overwrite_data=False)
#     for i in range(2,31):
#         if i == 21:
#             diag = DiagTable()
#             diag.add_file('atmos_6hrly', 6, 'hours', time_units='days')

#             #Write out diagnostics need for vertical interpolation post-processing
#             diag.add_field('dynamics', 'ps', time_avg=True)
#             diag.add_field('dynamics', 'bk')
#             diag.add_field('dynamics', 'pk')
#             diag.add_field('dynamics', 'zsurf')

#             #Tell model which diagnostics to write
#             diag.add_field('mixed_layer', 't_surf', time_avg=True)
#             diag.add_field('dynamics', 'temp', time_avg=True)
#             diag.add_field('atmosphere',   'temp_2m', time_avg=True)
#             diag.add_field('atmosphere',   'sphum_2m', time_avg=True)
#             diag.add_field('atmosphere',   'u_10m', time_avg=True)
#             diag.add_field('atmosphere',   'v_10m', time_avg=True) 

#             run_exp2.diag_table = diag 
#         run_exp2.run(i, num_cores=NCORES, overwrite_data=False)
