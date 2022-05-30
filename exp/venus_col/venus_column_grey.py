""" 
This script configures a column model that uses Isca's columnwise physics routines. 

Optical depths set to be those of Venus. 

Any questions to Neil Lewis:  
neil.lewis@physics.ox.ac.uk 
"""


import os

import numpy as np

from isca import ColumnCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE


# column model only uses 1 core
NCORES = 1

# compile code 
base_dir = os.path.dirname(os.path.realpath(__file__))
cb = ColumnCodeBase.from_directory(GFDL_BASE)
cb.compile() 

# create an Experiment object to handle the configuration of model parameters
exp = Experiment('venus_column5_conv_junktest', codebase=cb)

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 360, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('column', 'ps', time_avg=True)
diag.add_field('column', 'bk')
diag.add_field('column', 'pk')
diag.add_field('atmosphere', 't_surf', time_avg=True)
diag.add_field('mixed_layer', 'flux_lhe', time_avg=True)
diag.add_field('mixed_layer', 'flux_t', time_avg=True)
#diag.add_field('column', 'sphum', time_avg=True)
diag.add_field('column', 'ucomp', time_avg=True)
diag.add_field('column', 'vcomp', time_avg=True)
diag.add_field('column', 'temp', time_avg=True)
diag.add_field('two_stream', 'swdn_toa', time_avg=True)
diag.add_field('two_stream', 'olr', time_avg=True)
diag.add_field('two_stream', 'swdn_sfc', time_avg=True)
diag.add_field('two_stream', 'lwup_sfc', time_avg=True)
diag.add_field('two_stream', 'lwdn_sfc', time_avg=True)
diag.add_field('two_stream', 'net_lw_surf', time_avg=True)
diag.add_field('two_stream', 'tdt_rad', time_avg=True)
diag.add_field('two_stream', 'tdt_solar', time_avg=True)
diag.add_field('two_stream', 'flux_rad', time_avg=True)
diag.add_field('two_stream', 'flux_lw', time_avg=True)
diag.add_field('two_stream', 'flux_sw', time_avg=True)
diag.add_field('atmosphere', 'dt_ug_diffusion', time_avg=True)
diag.add_field('atmosphere', 'dt_vg_diffusion', time_avg=True)
diag.add_field('atmosphere', 'dt_tg_diffusion', time_avg=True)
exp.diag_table = diag

#Empty the run directory ready to run
exp.clear_rundir()


etas = np.array([2.265e-6, 2.925e-6, 4.500e-6, 7.950e-6, 1.525e-5, 
3.070e-5, 6.355e-5, 1.324e-4, 2.710e-4, 5.400e-4, 1.043e-3, 1.934e-3, 3.454e-3, 
5.929e-3, 9.814e-3, 1.569e-2, 2.430e-2, 3.657e-2, 5.361e-2, 7.672e-2, 1.074e-1, 
1.472e-1, 1.979e-1, 2.610e-1, 3.373e-1, 4.273e-1, 5.299e-1, 6.420e-1, 7.577e-1, 
8.679e-1, 9.602e-1])
central_half_etas = (etas[:-1] + etas[1:])/2.
half_etas = np.concatenate(([0.],central_half_etas,[1.])) 
half_etas = half_etas.tolist()

#Define values for the 'core' namelist
exp.namelist = namelist = Namelist({
    'main_nml':{
     'days'   : 360,
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,
     'dt_atmos':4800,
     'current_date' : [1,1,1,0,0,0],
     'calendar' : 'no_calendar'
    },
    
    'atmosphere_nml': {
        'idealized_moist_model': True
    },

    'column_nml': {
        'lon_max': 1, # number of columns in longitude, default begins at lon=0.0
        'lat_max': 1, # number of columns in latitude, precise 
                      # latitude can be set in column_grid_nml if only 1 lat used. 
        'num_levels': 31,  # number of levels 
        'initial_sphum': 0.0, 
        'valid_range_t': [50., 1500.], 
        'vert_coord_option': 'input', 
        'reference_sea_level_press': 9.2e6
    },

    'column_grid_nml': { 
        'lat_value': np.rad2deg(np.arcsin(1/np.sqrt(3))) # set latitude to that which causes insolation in frierson p2 radiation to be insolation / 4. 
        #'global_average': True # don't use this option at the moment
    },

    # set initial condition, NOTE: currently there is not an option to read in initial condition from a file. 
    'column_init_cond_nml': {
        'initial_temperature': 600., # initial atmospheric temperature 
        'surf_geopotential': 0.0, # applied to all columns 
        'surface_wind': 0.#5. # as described above 
    },

    'idealized_moist_phys_nml': {
        'do_damping': False, # no damping in column model, surface wind prescribed 
        'turb':True,        # DONT WANT TO USE THIS, BUT NOT DOING SO IS STOPPING MIXED LAYER FROM WORKING
        'mixed_layer_bc':True, # need surface, how is this trying to modify the wind field? ****
        'simple_surface':False,
        'do_simple': True, # simple RH calculation 
        'roughness_mom': 3.21e-05, # DONT WANT TO USE THIS, BUT NOT DOING SO IS STOPPING MIXED LAYER FROM WORKING
        'roughness_heat':3.21e-05,
        'roughness_moist':3.21e-05,                
        'two_stream_gray': True,     #Use grey radiation
        'convection_scheme': 'dry', #Use dry convection 
    },

    'two_stream_gray_rad_nml': {
        'rad_scheme': 'byrne',    
        'bog_a' :473.56 , # optical depth of background atmosphere 
        'bog_b': 0.0, # no condensible 
        'do_seasonal': False,                
        'atm_abs': 2.4,  
        'solar_constant': 2613.9 * (1 - 0.76), 
        'solar_exponent': 4.0, 
        'do_normal_integration_method': False 
    },
    
    
    'constants_nml': { # set constants to be those of Venus 
        'omega': 2.99e-7, 
        'grav': 8.87,
        'radius': 6.0518e6,
        'rdgas': 188.,
        'PSTD':92.0e6,
        'PSTD_MKS':92.0e5,
        'kappa': 188. / 850.1,
        #'CP_AIR': 900.,
    },

    
    'lscale_cond_nml': {
        'do_simple':True, # only rain 
        'do_evap':False,  # no re-evaporation of falling precipitation 
    },

    'surface_flux_nml': {
        'use_virtual_temp': False, # use virtual temperature for BL stability 
        'do_simple': True,
        'old_dtaudv': True    
    },

    'vert_turb_driver_nml': { # DONT WANT TO USE THIS, BUT NOT DOING SO IS STOPPING MIXED LAYER FROM WORKING
        'do_mellor_yamada': False,     # default: True
        'do_diffusivity': False,        # default: False
        'do_simple': True,             # default: False
        'constant_gust': 0.0,          # default: 1.0
        'use_tau': False
    },

    #Use a large mixed-layer depth, and the Albedo of the CTRL case in Jucker & Gerber, 2017
    'mixed_layer_nml': {
        'tconst' : 600.,
        'prescribe_initial_dist':False,
        'evaporation':False,   
        'depth': 5.,                          #Depth of mixed layer used
        'albedo_value': 0.0,                  #Albedo value used             
    },

    'dry_convection_nml': {
        'tau':7200.,
        'Gamma':1.0, 

    },
    
    'sat_vapor_pres_nml': {
        'do_simple':True, 
        'tcmin_simple' : -200, 
        'tcmax_simple' : 2000
    },

    # define pressure coordinate 
    'vert_coordinate_nml': {
        'bk': half_etas,
        'pk': np.zeros_like(half_etas).tolist(),
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
    exp.run(1, use_restart=False, num_cores=NCORES, overwrite_data=False)
    for i in range(2,151):
        exp.run(i, num_cores=NCORES, overwrite_data=False)


