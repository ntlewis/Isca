"""
GROUP1: 
om = 1/4,  a = 2-1/8
om = 1/16, a = 2-1/8 
+ 
om = 1,    a = 1/8
om = 1/2,  a = 1/4
om = 1/8,  a = 1 
+ 
om = 1/8,  a = 1/4 
om = 1/32, a = 1 
om = 1/64, a = 2

with: ps = 1bar and no sw absorption 
"""

import os

import numpy as np

from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE


#### set constants and experiment arrays 

# constants 
pse            = 100000.                                  # surface pressure
grav           = 9.81                                     # acceleration due to gravity 
tau_lw_s       = 2.                                       # longwave optical depth 
tau_sw_s       = 2./3. * 0.25 * 160. * tau_lw_s *12./16.  # shortwave optical depth in non-window region 
solar_constant = 1360.                                    # solar constant 
albedo         = 0.27                                     # bond albedo 
omegae         = 7.29e-5                                  # Earth's rotation rate 
radiuse        = 6.371e6                                  # Earth's radius (m)
R_d            = 287.                                     # gas constant for dry air 
cp             = 1005.                                    # specific heat capacity at constant pressure 
dt_conv        = 3600. * 2.                               # timescale for dry convection 
gamma_conv     = 0.7                                      # sets lapse rate in dry convection scheme: adjust to stable profile (check sensitivity to this)
ml_depth       = 1.                                       # mixed layer depth (m)
pr             = 20000.                                   # reference pressure for optical depths (unused)
lw_f           = 0.1                                      # longwave p-dependence transition fraction (unused)
sw_window_frac = 1.                                       # shortwave window fraction (if 1, atmosphere is transparent)
sw_window_flag = True                                     # flag for shortwave window 

# experiment arrays 
om_array  = np.array([1./4.])#, 1./4., 1./4., 1./4., 1./4., 1./16., 1./16., 1./16., 1./16., 1./16., 1.   , 1./2., 1./8., 1./8., 1./32., 1./64])
a_array   = np.array([2.   ])#, 1.   , 1./2., 1./4., 1./8., 2.    , 1.    , 1./2. , 1./4. , 1./8. , 1./8., 1./4., 1.   , 1./4., 1.    , 2.   ])
ps_array  = np.array([1.   ])#, 1.   , 1.   , 1.   , 1.   , 1.    , 1.    , 1.    , 1.    , 1.    , 1.   , 1.   , 1.   , 1.   , 1.    , 1.   ])
abs_array = np.array([0.   ])#, 0.   , 0.   , 0.   , 0.   , 0.    , 0.    , 0.    , 0.    , 0.    , 0.   , 0.   , 0.   , 0.   , 0.    , 0.   ])

# damping timescale depends on radius 
tsc_array = np.array([16.4 ])#, 12.2 , 9.1  , 6.8  , 5.1  , 16.4  , 12.2  , 9.1   , 6.8   , 5.1   , 5.1  , 6.8  , 12.2 , 6.8  , 12.2  , 16.4 ])

# names corresponding to experiment arrays 
om_names  =          ['1_4']#, '1_4', '1_4', '1_4', '1_4', '1_16', '1_16', '1_16', '1_16', '1_16', '1_1', '1_2', '1_8', '1_8', '1_32', '1_32']
a_names   =          ['2_1']#, '1_1', '1_2', '1_4', '1_8', '2_1' , '1_1' , '1_2' , '1_4' , '1_8' , '1_8', '1_4', '1_1', '1_4', '1_1' , '1_2' ]
ps_names  =          ['1'  ]#, '1'  , '1'  , '1'  , '1'  , '1'   , '1'   , '1'   , '1'   , '1'   , '1'  , '1'  , '1'  , '1'  , '1'   , '1'   ]
abs_names =          ['0'  ]#, '0'  , '0'  , '0'  , '0'  , '0'   , '0'   , '0'   , '0'   , '0'   , '0'  , '0'  , '0'  , '0'  , '0'   , '0'   ] 


###### SETTING UP THE MODEL #######
# set up base experiment including diagnostics and namelist 

NCORES = 32
base_dir = os.path.dirname(os.path.realpath(__file__))
cb = IscaCodeBase.from_directory(GFDL_BASE)
cb.compile()  

# template file name and experiment 
exp = Experiment('grey_T42L50_om_X_Y_a_X_Y_ps_Zbar_abs_X_Yfrac', codebase=cb)


# format for diagnostics 
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')

# which diagnostics to write (will get more out / higher frequency output for experiments of interest using future runs)
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'omega', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'height', time_avg=True)
diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('two_stream_scatter', 'olr', time_avg=True)
diag.add_field('two_stream_scatter', 'swdn_toa', time_avg=True)
#diag.add_field('damping', 'udt_rdamp', time_avg=True)

exp.diag_table = diag

#Empty the run directory ready to run
exp.clear_rundir()

#Define values for the 'core' namelist
exp.namelist = Namelist({
    'main_nml':{
     'days'   : 360,
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,
     'dt_atmos':30,
     'current_date' : [1,1,1,0,0,0],
     'calendar' : 'thirty_day'
    },

    'atmosphere_nml': {
        'idealized_moist_model': True
    },

    'spectral_dynamics_nml': {
        'damping_order'            : 2,                 
        'damping_coeff'            : 1. / (12.2 * 60. * 60.), 
        # laplacian sponge 
        'eddy_sponge_coeff'        : 1. / (0.4 * 60. * 60.),
        'water_correction_limit'   : 200.e2,
        'reference_sea_level_press': pse,
        'valid_range_t'            : [100.,800.],
        'initial_sphum'            : [0.0],
        'robert_coeff'             : 0.03,
        # T42 RESOLUTION
        'lon_max'                  : 128,
        'lat_max'                  : 64,
        'num_fourier'              : 42,
        'num_spherical'            : 43,
        # 50 vertical levels with high top and sufficient resolution in stratosphere 
        'vert_coord_option'        : 'hybrid',        
        'num_levels'               : 50,
        'exponent'                 : 7.0, 
        'surf_res'                 : 0.5, 
        'scale_heights'            : 11.0, 
        # INITIAL CONDITION to be taken from column model. File specified in experiment specific loop. 
        'initial_state_option'    :'input' # quiescent for cold start 
    },


    'ic_from_external_file_nml':{  # file for initial condition specified in experiment specific loop.
        'u_name': 'ucomp', 
        'v_name': 'vcomp', 
        't_name': 'temp', 
        'ps_name': 'ps'
    },
    

    'idealized_moist_phys_nml': {
        # no damping: sponge is laplacian in dynamical core 
        'do_damping'        : False, 
        'turb'              : True,
        'mixed_layer_bc'    : True,
        'do_virtual'        : False,
        'do_simple'         : True,
        'roughness_mom'     : 5.e-3,
        'roughness_heat'    : 5.e-3,
        'roughness_moist'   : 1.e-5,                
        'two_stream_gray'   : False,   
        # use my two_stream_scatter module for rt 
        'two_stream_scatter': True, 
        # use dry convection scheme 
        'convection_scheme' : 'dry', 
    },

  

    'two_stream_scatter_rad_nml': {
        'sw_optical_depth': 'generic', 
        'lw_abs_a'        : 0.0, 
        'lw_abs_b'        : 0.0, 
        'lw_abs_c'        : 0.0, 
        'lw_abs_d'        : 0.0,
        # longwave optical depth 
        'lw_abs_e'        : tau_lw_s, 
        'lw_abs_f'        : 0.0, 
        'lw_abs_pref'     : pr, #(unused)
        'lw_sca_a'        : 0.0,
        'lw_sca_b'        : 0.0,
        'sw_abs_a'        : 0.0,
        'sw_abs_b'        : 0.0, 
        # shortwave optical depth in non-window region 
        'sw_abs_c'        : tau_sw_s, 
        'sw_sca_a'        : 0.0, 
        'sw_sca_b'        : 0.0,
        # shortwave window fraction 
        'sw_window_frac'  : sw_window_frac, 
        'do_sw_window'    : sw_window_flag, 
        # assymetry factor for scattering (unused if sw_sca's are 0.0)
        'g_asym'          : 0.0, 
        # do true perpetual equinox, not frierson p2 
        'do_seasonal'     : False,
        'perp_eq'         : True, 
        # apply albedo here, otherwise things get complicated with the window
        'solar_constant'  : solar_constant * (1 - albedo)
    },

    'vert_turb_driver_nml': {
        'do_mellor_yamada': False,    
        'do_diffusivity'  : True,       
        'do_simple'       : True,          
        'constant_gust'   : 1.0,         
        'use_tau'         : False
    },
    
    'diffusivity_nml': {
        'do_entrain': False,
        'do_simple' : True,
    },

    'surface_flux_nml': {
        'use_virtual_temp': False,
        'do_simple'       : True,
        'old_dtaudv'      : True    
    },

    'constants_nml':{
        # set planet constants 
        'omega'   : omegae, # modified in experiment loop 
        'grav'    : grav, 
        'radius'  : radiuse, # modified in experiment loop 
        'rdgas'   : R_d, 
        'PSTD'    : pse * 10., 
        'PSTD_MKS': pse, 
        'kappa'   : R_d / cp,
    },



    'mixed_layer_nml': {
        # use initial condition (file provided in experiment loop)
        'ic_neil_style': True, 
        # no evap because simulation is dry 
        'evaporation'  : False,   
        'depth'        : ml_depth, 
        # albedo zero because applied to solar constant instead                        
        'albedo_value' : 0.0,                          
    },

    'dry_convection_nml': {
        'tau'  : dt_conv,
        'Gamma': gamma_conv, 

    },
    
    'lscale_cond_nml': {
        'do_simple':True,
        # no evaporation because experiments are dry 
        'do_evap': False 
    },
    
    'sat_vapor_pres_nml': {
        'do_simple'    : True,
        # range here doesn't matter as experiments are dry, but this prevents the model from crashing 
        'tcmin_simple' : -273, 
        'tcmax_simple' : 2000
    },




    # FMS Framework configuration
    'diag_manager_nml': {
        'mix_snapshot_average_fields': False  
    },

    'fms_nml': {
        'domains_stack_size': 600000                    
    },

    'fms_io_nml': {
        'threading_write': 'single',                       
        'fileset_write'  : 'single',                  
    },


})

# Run experiments 
if __name__=="__main__":


    # loop over experiments 
    for n, om in enumerate(om_array): 
        # derive experiment 
        run_exp = exp.derive('grey_T42L50_om_'+om_names[n]+'_a_'+a_names[n]+'_ps_'+ps_names[n]+'bar_abs_'+abs_names[n]+'frac_dt30')

        # modify namelist 

        # constants 
        run_exp.namelist['constants_nml']['radius']   = radiuse * a_array[n] 
        run_exp.namelist['constants_nml']['omega']    = omegae * om_array[n] 
        run_exp.namelist['constants_nml']['PSTD']     = pse * ps_array[n] * 10. 
        run_exp.namelist['constants_nml']['PSTD_MKS'] = pse * ps_array[n] 
        run_exp.namelist['spectral_dynamics_nml']['reference_sea_level_press'] = pse * ps_array[n]
        run_exp.namelist['two_stream_scatter_rad_nml']['sw_window_frac'] = 1. - abs_array[n]

        # hyperdiffusion timescale (to scale with radius)
        run_exp.namelist['spectral_dynamics_nml']['damping_coeff'] = 1. / (tsc_array[n] * 60. * 60.)

        # initial condition related things 
        run_exp.inputfiles = [os.path.join(base_dir,'/home/lewis/research/isca_dir/Isca/exp/grey_suprot/column_ics/ic_files/'+ps_names[n]+'bar_'+abs_names[n]+'frac_ic.nc')]
        run_exp.namelist['ic_from_external_file_nml']['file_name'] = 'INPUT/'+ps_names[n]+'bar_'+abs_names[n]+'frac_ic.nc'
        run_exp.namelist['mixed_layer_nml']['neil_ic_file_name']   = 'INPUT/'+ps_names[n]+'bar_'+abs_names[n]+'frac_ic.nc'


        # actually run experiments 
        run_exp.run(1, use_restart=False, num_cores=NCORES, overwrite_data=False)
        for i in range(2,5):
            run_exp.run(i, num_cores=NCORES, overwrite_data=False)
