import os

import numpy as np

from isca import SocratesCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE, FailedRunError

NCORES = 16
base_dir = os.path.dirname(os.path.realpath(__file__))
# a CodeBase can be a directory on the computer,
# useful for iterative development
cb = SocratesCodeBase.from_directory(GFDL_BASE)

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
exp = Experiment('nakajima_eccentricity', codebase=cb)
#exp = Experiment('byrne_sc_test', codebase=cb)



#Empty the run directory ready to run
exp.clear_rundir()

exp.inputfiles=['sphum_input.nc']

#Define values for the 'core' namelist
exp.namelist = namelist = Namelist({
    'main_nml':{
     'days'   : 360,
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,
     'dt_atmos':600.,
     'current_date' : [0001,1,1,0,0,0],
     'calendar' : 'thirty_day'
    },

    'idealized_moist_phys_nml': {
        'do_damping': True,
        'turb':True,
        'mixed_layer_bc':True,
        'do_virtual' :False,
        'do_simple': False,
        'roughness_mom':3.21e-05,
        'roughness_heat':3.21e-05,
        'roughness_moist':3.21e-05,                
        'two_stream_gray': False,     #Use grey radiation
	'do_rrtm_radiation': False, 
	'do_socrates_radiation':True, 
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

    'surface_flux_nml': {
        'use_virtual_temp': False,
        'do_simple': False,
        'old_dtaudv': True    
    },

    'atmosphere_nml': {
        'idealized_moist_model': True
    },

    #Use a large mixed-layer depth, and the Albedo of the CTRL case in Jucker & Gerber, 2017
    'mixed_layer_nml': {
        'tconst' : 290.,
        'prescribe_initial_dist':True,
        'evaporation':False,   
        'depth': 20.0,                          #Depth of mixed layer used
        'albedo_value': 0.2,                  #Albedo value used             
        'do_qflux': False,
    },

    'dry_convection_nml': {
        'tau':7200.,
        'Gamma':1.0, 

    },

    'betts_miller_nml': {
       'rhbm': .7   , 
       'do_simp': False, 
       'do_shallower': True, 
    },
    
    'lscale_cond_nml': {
        'do_simple':False,
        'do_evap':True
    },
    
    'sat_vapor_pres_nml': {
        'do_simple':False,
	'tcmin': -240., 
	'tcmax':  200.,
	
    },
    
    'damping_driver_nml': {
        'do_rayleigh': True,
        'trayfric': -0.25,              # neg. value: time in *days*
        'sponge_pbottom':  5000.,           #Bottom of the model's sponge down to 50hPa (units are Pa)
        'do_conserve_energy': True,             
    },

    'socrates_rad_nml': {
        'stellar_constant':1370.,
        'lw_spectral_filename': '/scratch/nl290/ROCKE/sp_lw_ga7/sp_lw_ga7_dsa',
        'sw_spectral_filename': '/scratch/nl290/ROCKE/sp_sw_ga7/sp_sw_ga7_dsa_sun',
        'tidally_locked': False,
        'do_read_ozone': False,
        'ozone_file_name':'ozone_1990',
        'ozone_field_name':'ozone_1990',
        'do_read_h2o': True, 
        'h2o_file_name': 'sphum_input', 
        'h2o_field_name': 'sphum_input', 
        'dt_rad': 7200,
        'store_intermediate_rad': True,
        'chunk_size': 16,
        'use_pressure_interp_for_half_levels':False,
	'equinox_day':0.0,
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
        'num_levels':40,               #How many model pressure levels to use
        'valid_range_t':[40.,800.],
        'initial_sphum':[0.0],
        'vert_coord_option':'uneven_sigma', #Use the vertical levels from Frierson 2006
        'surf_res':0.2,
        'scale_heights' : 11.0,
        'exponent':7.0,
        'robert_coeff':0.03, 
        'do_water_correction': False 
    },

    'astronomy_nml': {
        'ecc': 0.0,
        'obliq': 0.0,     # planetary obliquity
        'per': 0.0        # begin at autumn equinox
    },
})

#Lets do a run!
#eccentricities = np.array([0.8])
eccentricities = np.array([0.0, 0.15, 0.3, 0.45, 0.6])#, 0.75])#, 0.75, 0.9])#0.6, 0.75, 0.9])#np.array([0.8,0.6,0.4,0.2,0.0])
if __name__=="__main__":
        resolution = 'T21', 40
        exp.set_resolution(*resolution)
	for e in eccentricities:
		run_exp = exp.derive('soc_eccentricity_'+str(e)+'_20_new_dry')
		run_exp.namelist['astronomy_nml']['ecc'] = e
		# also change semi-maj axis
		run_exp.namelist['socrates_rad_nml']['stellar_constant'] = 1368. * (1 - e**2.)**(1./2.)

                #Tell model how to write diagnostics
                diag = DiagTable()
                diag.add_file('atmos_monthly', 5, 'days', time_units='days')
                #Tell model which diagnostics to write
                diag.add_field('dynamics', 'ps', time_avg=True)
                diag.add_field('atmosphere', 'rh', time_avg=True)
                diag.add_field('dynamics', 'bk')
                diag.add_field('dynamics', 'pk')
                diag.add_field('atmosphere', 'precipitation', time_avg=True)
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
		diag.add_field('socrates', 'soc_olr', time_avg=True)
		diag.add_field('socrates', 'soc_toa_sw', time_avg=True)
		diag.add_field('socrates', 'soc_surf_flux_sw', time_avg=True)
		diag.add_field('socrates', 'soc_surf_flux_lw', time_avg=True)
                run_exp.diag_table = diag
                
		run_exp.run(1, use_restart=False, num_cores=NCORES, overwrite_data=False)
                
                
                
		for i in range(2,21):
                        if i == 9:
                                diag = DiagTable()
                                diag.add_file('atmos_monthly', 2, 'days', time_units='days')
                                #Tell model which diagnostics to write
                                diag.add_field('dynamics', 'ps', time_avg=True)
                                diag.add_field('atmosphere', 'rh', time_avg=True)
                                diag.add_field('dynamics', 'bk')
                                diag.add_field('dynamics', 'pk')
                                diag.add_field('atmosphere', 'precipitation', time_avg=True)
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
				diag.add_field('socrates', 'soc_olr', time_avg=True)
				diag.add_field('socrates', 'soc_toa_sw', time_avg=True)
				diag.add_field('socrates', 'soc_surf_flux_sw', time_avg=True)
				diag.add_field('socrates', 'soc_surf_flux_lw', time_avg=True)
                                run_exp.diag_table = diag
                                
			run_exp.run(i, num_cores=NCORES, overwrite_data=False)

    #for i in range(2,121):
    #    exp.run(i, num_cores=NCORES)
