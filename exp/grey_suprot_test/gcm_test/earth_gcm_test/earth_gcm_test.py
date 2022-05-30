""" 
Simple grey titan model 
"""


import os

import numpy as np

from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE


# column model only uses 1 core
NCORES = 32

# compile code 
base_dir = os.path.dirname(os.path.realpath(__file__))
cb = IscaCodeBase.from_directory(GFDL_BASE)
cb.compile() 

# create an Experiment object to handle the configuration of model parameters
exp = Experiment('earth_gcm_test', codebase=cb)

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('mixed_layer', 'flux_t', time_avg=True)
#diag.add_field('column', 'sphum', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('two_stream_scatter', 'swdn_toa', time_avg=True)
diag.add_field('two_stream_scatter', 'olr', time_avg=True)
exp.diag_table = diag

#Empty the run directory ready to run
exp.clear_rundir()

#list planets parameters here 
ps = 100000
T_ini = 270.
grav = 9.81
tau_lw_s = 3./2. * 3. * 1.25 #was4 #(corresponds to 50,000 at surface) #473.56
tau_sw_s = 2./3. * 0.25 * 160. * tau_lw_s * 12./16.  #.5#0.22 * tau_lw_s #0.22 is anti-greenhouse parameter 
solar_constant = 1360.
albedo = 0.27
omega = 7.29e-5
radius = 6.371e6
R_d = 287.
cp = 1005.
dt_conv = 3600. * 4. # slow timescale for dry convection 
ml_depth = 5. # meters 
sw_window_frac = 1. #- 4./7.
sw_window_flag = True 






#Define values for the 'core' namelist
exp.namelist = namelist = Namelist({
    'main_nml':{
     'days'   : 360,
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,
     'dt_atmos':720,
     'current_date' : [1,1,1,0,0,0],
     'calendar' : 'thirty_day'
    },
    
    'atmosphere_nml': {
        'idealized_moist_model': True
    },

    
    'spectral_init_cond_nml': {
        'initial_temperature': T_ini, # initial atmospheric temperature 
    },

    'spectral_dynamics_nml': {
        'damping_order': 4,       
        'damping_coeff':1 / (1./6. * 86400.),
        'water_correction_limit': 0.0,
        'reference_sea_level_press':ps,
        'num_levels':45,               #How many model pressure levels to use
        'valid_range_t':[10.,800.],
        'initial_sphum':[0.0],
        'vert_coord_option':'uneven_sigma', #Use the vertical levels from Frierson 2006
        'surf_res':0.1,
        'scale_heights' :17,
        'exponent':2.5,
        'robert_coeff':0.03, #        'raw_filter_coeff': 0.0,#53, 
        'lon_max': 128,
        'lat_max': 64,
        'num_fourier': 42,
        'num_spherical': 43,

    },

    'idealized_moist_phys_nml': {
        'do_damping': True, 
        'turb':True,        
        'mixed_layer_bc':True, 
        'do_simple': True, # simple RH calculation 
        'roughness_mom': 3.21e-05, # DONT WANT TO USE THIS, BUT NOT DOING SO IS STOPPING MIXED LAYER FROM WORKING
        'roughness_heat':3.21e-05,
        'roughness_moist':3.21e-05,                
        'two_stream_gray': False,     #Use grey radiation
        'two_stream_scatter':True, 
        'do_rrtm_radiation': False, 
        'do_socrates_radiation': False, 
        'convection_scheme': 'dry', #Use the simple Betts Miller convection scheme 
    },

    'two_stream_scatter_rad_nml': {
        'sw_optical_depth': 'generic', 
        'lw_abs_a': 0.0, 
        'lw_abs_b': 0.0, 
        'lw_abs_c': 0.0, 
        'lw_abs_d': 0.0,#tau_lw_s, #grav * tau_lw_s / ps, 
        'lw_abs_e': grav * tau_lw_s / ps, 
        'lw_sca_a':0.0,
        'lw_sca_b':0.0,
        'sw_abs_a': 0.0,#grav * tau_sw_s / ps, 
        'sw_abs_b':0.0, 
        'sw_abs_c':grav * tau_sw_s / ps, 
        'sw_sca_a':0.0, 
        'sw_sca_b':0.0,
        'sw_window_frac':sw_window_frac, 
        'do_sw_window':sw_window_flag, 
        'g_asym' : 0.0, 
        'do_seasonal':False,
        'solar_constant': solar_constant * (1 - albedo) #2613.9 * (1 - 0.76)
    },

    'constants_nml': { # set constants to be those of Titan
        'omega': omega, 
        'grav': grav, 
        'radius': radius, 
        'rdgas': R_d, 
        'PSTD': ps * 10., #dynes cm-2
        'PSTD_MKS': ps, #pa 
        'kappa': R_d / cp,
        #'CP_AIR': 900.,
    },

    'astronomy_nml':{
        'obliq':0.0, 
        'ecc':0.0,
    },


    # 'llcs_nml': { 
    #     'llcs_rhcrit': 0.8, 
    #     'llcs_timescale': dt_conv, 
    #     'cloud_option': 'ALL_RAIN', 
    #     'do_simple': True, 
    # },

    'dry_convection_nml': {
        'tau':dt_conv,
        'Gamma':1.0, 

    },
    
    'lscale_cond_nml': {
        'do_simple':True, # only rain 
        'do_evap':False,  # no re-evaporation of falling precipitation 
    },
    
    'damping_driver_nml': {
        'do_rayleigh': True,
        'rayleigh_eddy_only':False, 
        'trayfric': -0.3,              # neg. value: time in *days*
        'sponge_pbottom':  0.5 * (ps/1.5e5),           #Bottom of the model's sponge down to 50hPa (units are Pa)
        'do_conserve_energy': False,             
    },


    'surface_flux_nml': {
        'use_virtual_temp': False, # use virtual temperature for BL stability 
        'do_simple': True,
        'old_dtaudv': True    
    },

    'vert_turb_driver_nml': { # DONT WANT TO USE THIS, BUT NOT DOING SO IS STOPPING MIXED LAYER FROM WORKING
        'do_mellor_yamada': False,     # default: True
        'do_diffusivity': True,        # default: False
        'do_simple': True,             # default: False
        'constant_gust': 0.1,          # default: 1.0
        'use_tau': False
    },

    #Use a large mixed-layer depth, and the Albedo of the CTRL case in Jucker & Gerber, 2017
    'mixed_layer_nml': {
        'tconst' : T_ini,
        'prescribe_initial_dist':False,
        'evaporation':False,   
        'depth': ml_depth,                          #Depth of mixed layer used
        'albedo_value': 0.0,                  #Albedo value used             
    },
    
    'sat_vapor_pres_nml': {
        'do_simple':True, 
        'tcmin_simple' : -273, 
        'tcmax_simple' : 2000
    },

    # # define pressure coordinate 
    # 'vert_coordinate_nml': {
    #     'bk': [0.000000, 0.0117665, 0.0196679, 0.0315244, 0.0485411, 0.0719344, 0.1027829, 0.1418581, 0.1894648, 0.2453219, 0.3085103, 0.3775033, 0.4502789, 0.5244989, 0.5977253, 0.6676441, 0.7322627, 0.7900587, 0.8400683, 0.8819111, 0.9157609, 0.9422770, 0.9625127, 0.9778177, 0.9897489, 1.0000000],
    #     'pk': [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000],
    #    },

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
    for i in range(2,11):
        exp.run(i, num_cores=NCORES, overwrite_data=False)
