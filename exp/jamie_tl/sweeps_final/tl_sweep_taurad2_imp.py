import numpy as np

from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

cb = IscaCodeBase.from_directory(GFDL_BASE)

#expname
expname = 'mv_T42L20_om_1_8_a_1_1_ps_1_1_spg'
timestep  = 600#60#60.#120

# planet params 
omega          = 7.29e-5 
radius         = 6.371e6 #/8.#/ 20.
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
spongep = 0.05*ps

# ml params 
mldepth = 2.5 
mlalbedo = 0.0


# core params 
fourier   = 42 
spherical = fourier+1
lats      = 64
lons      = lats * 2
cores     = 32
levels    = 25#48
T_ini     = 280. 

cb.compile()  

exp = Experiment(expname, codebase=cb)

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 100., 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'height', time_avg=True)
diag.add_field('dynamics', 'omega', time_avg=True)
diag.add_field('hs_forcing', 'teq', time_avg=True)
#diag.add_field('atmosphere', 'dt_tg_diffusion', time_avg=True)
#diag.add_field('atmosphere', 'dt_tg_convection', time_avg=True)

exp.diag_table = diag

#Empty the run directory ready to run
exp.clear_rundir()

# define namelist values as python dictionary
# wrapped as a namelist object.
namelist = Namelist({
    'main_nml': {
     'days'   : 100,
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,
     'dt_atmos':timestep,
     'current_date' : [1,1,1,0,0,0],
     'calendar' : 'no_calendar'
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
        'no_diffm_tendency':True, 
        'no_difft_tendency':True
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
    'hs_forcing_nml': {
        'equilibrium_t_option':'EXOPLANET_TL',
        't_zero': 365.,    # temperature at reference pressure at equator (default 315K)
        't_strat': 220.,   # stratosphere temperature (default 200K)
        'delh': 145.,       # equator-pole temp gradient (default 60K)
        'delv': 10.,       # lapse rate (default 10K)
        'eps': 0.,         # stratospheric latitudinal variation (default 0K)
        'sigma_b': 0.7,    # boundary layer friction height (default p/ps = sigma = 0.7)

        # negative sign is a flag indicating that the units are days
        'ka':   -6.,#*32*np.sqrt(2),      # Constant Newtonian cooling timescale (default 40 days)
        'ks':    -1.5,#*32*np.sqrt(2),      # Boundary layer dependent cooling timescale (default 4 days)
        'kf':   -1./1.,       # BL momentum frictional timescale (default 1 days)

        'do_conserve_energy':   True,  # convert dissipated momentum into heat (default True)
    },

    
    
    'spectral_dynamics_nml': {
        'damping_order'           : 4,                      # default: 2
        'water_correction_limit'  : 200.e2,                 # default: 0
        'reference_sea_level_press': 1.0e5,                  # default: 101325
        'valid_range_t'           : [100., 800.],           # default: (100, 500)
        'initial_sphum'           : 0.0,                  # default: 0
        'damping_coeff':10./86400.*1.,
        #'eddy_sponge_coeff':10./86400. * radius**2. / (42*43),
        #'zmu_sponge_coeff': 10./86400. * radius**2. / (42*43),
        #'zmv_sponge_coeff': 10./86400. * radius**2. / (42*43),
        'vert_coord_option'       : 'hybrid',         # default: 'even_sigma'
        'scale_heights': 5.0,
        'exponent': 3.0,
        'surf_res': 0.05,
        'num_levels':25, 
        'lat_max': 64,
        'lon_max': 128,
        'num_fourier':42,
        'num_spherical':42+1
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

om = 16. 
rads = np.array([4.,2.,1.,1./2.,1./4.])#1. / 1.
names = ['4_1tau','2_1tau','1_1tau','1_2tau','1_4tau']
radnames = ['4_1rad', '2_1rad', '1_1rad', '1_2rad', '1_4rad']
taus  = np.array([4.,2.,1.,0.5,0.25])
#steps = np.array([64.,64.,64.,64.])#32.,32.,32.,64.,128.,256.,512.,512.])
steps = np.array([600.,600.,600.,300.,180.])
#Lets do a run!
if __name__ == '__main__':

    for n, tau in enumerate(taus):
        for m, rad in enumerate(rads):
            if n == 2 and m == 2:
                continue
            omega_exp = omega / om
            rad_exp = radius * rad
            ka_exp = -6. * tau 
            ks_exp = -1.5 * tau 
            step_exp = steps[m]
            run_exp = exp.derive('F_tl_'+radnames[m]+'_'+'1_16om_'+names[n]+'_rad')
            run_exp.namelist['main_nml']['dt_atmos'] = step_exp
            run_exp.namelist['constants_nml']['omega'] = omega_exp
            run_exp.namelist['constants_nml']['radius'] = rad_exp 
            run_exp.namelist['hs_forcing_nml']['ka'] = ka_exp 
            run_exp.namelist['hs_forcing_nml']['ks'] = ks_exp
            run_exp.run(1, num_cores=cores,use_restart=False,overwrite_data=True)
            for i in range(2, 31):
                if i == 28: 
                    diag = DiagTable()
                    diag.add_file('atmos_monthly', 1., 'days', time_units='days')
                    #Tell model which diagnostics to write
                    diag.add_field('dynamics', 'ps', time_avg=True)
                    diag.add_field('dynamics', 'bk')
                    diag.add_field('dynamics', 'pk')
                    diag.add_field('dynamics', 'ucomp', time_avg=True)
                    diag.add_field('dynamics', 'vcomp', time_avg=True)
                    diag.add_field('dynamics', 'temp', time_avg=True)
                    diag.add_field('dynamics', 'height', time_avg=True)
                    diag.add_field('dynamics', 'omega', time_avg=True)
                    diag.add_field('hs_forcing', 'teq', time_avg=True)
                    run_exp.diag_table = diag
                run_exp.run(i, num_cores=cores, overwrite_data=True) 