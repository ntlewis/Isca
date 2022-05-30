import os

import numpy as np

from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES = 24
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
exp = Experiment('fr_T31_L20_dt240_even_rad_1_8_damp_1_10', codebase=cb)

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('atmosphere', 'precipitation', time_avg=True)
diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('dynamics', 'sphum', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'vor', time_avg=True)
diag.add_field('dynamics', 'div', time_avg=True)
#diag.add_field('damping', 'udt_rdamp', time_avg=True)

exp.diag_table = diag

#Empty the run directory ready to run
exp.clear_rundir()

#Define values for the 'core' namelist
exp.namelist = namelist = Namelist({
    'main_nml':{
     'days'   : 30,
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,
     'dt_atmos':240,
     'current_date' : [1,1,1,0,0,0],
     'calendar' : 'thirty_day'
    },

    'atmosphere_nml': {
        'idealized_moist_model': True
    },

    'spectral_dynamics_nml': {
        'damping_order': 4,                   
        'damping_coeff': 1 / (1./1. * 86400.), #'eddy_sponge_coeff' : 1/(1./10. * 86400.),
        'water_correction_limit': 200.e2,
        'reference_sea_level_press':1.0e5,
        'num_levels':20,#40,#200,               #How many model pressure levels to use
        'valid_range_t':[100.,800.],
        'initial_sphum':[0.0],
        'vert_coord_option':'even_sigma',#'hybrid', #Use the vertical levels from Frierson 2006
        'robert_coeff':0.03,
        'lon_max': 96,
        'lat_max': 48,
        'num_fourier': 31,
        'num_spherical': 32,
        #'alpha_implicit':1.0,
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
        'two_stream_gray': True,     #Use grey radiation
        'convection_scheme': 'dry', #Use the simple Betts Miller convection scheme from Frierson
    },

    'vert_turb_driver_nml': {
        'do_mellor_yamada': False,     # default: True
        'do_diffusivity': True,        # default: False
        'do_simple': True,             # default: False
        'constant_gust': 1.0,          # default: 1.0
        'use_tau': False
    },
    
    'diffusivity_nml': {
        'do_entrain':False,
        'do_simple': True,
        #'background_m':0.1, 
        #'background_t':0.1,
    },

    'surface_flux_nml': {
        'use_virtual_temp': False,
        'do_simple': True,
        'old_dtaudv': True    
    },

    'constants_nml':{
        'radius': 0.125 * 6400.0e3
    },



    #Use a large mixed-layer depth, and the Albedo of the CTRL case in Jucker & Gerber, 2017
    'mixed_layer_nml': {
        'tconst' : 285.,
        'prescribe_initial_dist':True,
        'evaporation':False,   
        'depth': 1.,                          #Depth of mixed layer used
        'albedo_value': 0.31,                  #Albedo value used             
    },

    'dry_convection_nml': {
        'tau':2.*3600.,
        'Gamma':0.7, 

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

    
    'damping_driver_nml': {
        'do_rayleigh': False,
        'trayfric': -(1./2.),              # neg. value: time in *days*
        'sponge_pbottom':  50.,           #Bottom of the model's sponge down to 50Pa (units are Pa)
        'do_conserve_energy': False,            
        'rayleigh_eddy_only':False, 
    },

    'two_stream_gray_rad_nml': {
        'rad_scheme': 'frierson',            #Select radiation scheme to use, which in this case is Frierson
        'do_seasonal': False,                #do_seasonal=false uses the p2 insolation profile from Frierson 2006. do_seasonal=True uses the GFDL astronomy module to calculate seasonally-varying insolation.
        'atm_abs': 0.0,                      # default: 0.0        
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


})

#Lets do a run!
if __name__=="__main__":
    exp.run(1, use_restart=False, num_cores=NCORES, overwrite_data=True)
    for i in range(2,37):
        exp.run(i, num_cores=NCORES, overwrite_data=True)
