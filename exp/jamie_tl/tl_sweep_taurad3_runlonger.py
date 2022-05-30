import numpy as np

from isca import DryCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES = 32
#RESOLUTION = 'T42', 20  # T42 horizontal resolution, 25 levels in pressure

omega_e = 7.29e-5 
omega  = omega_e /16.#/16.#(32.*np.sqrt(2))#4. #$#/ 10. 
radius = 6.371e6#/8. #*2. #/8.
grav   = 9.81 
ps     = 1.e5 
Rd     = 287. 
cp     = 1005. 

# a CodeBase can be a directory on the computer,
# useful for iterative development
cb = DryCodeBase.from_directory(GFDL_BASE)

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

exp_name = 'tl_jamie_test_1rad_1_16om_25_dt300_tauchange1_diff4'#_T170L12'
exp = Experiment(exp_name, codebase=cb)

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
#diag.add_field('hs_forcing', 'teq', time_avg=True)
#diag.add_field('dynamics', 'vor', time_avg=True)
#diag.add_field('dynamics', 'div', time_avg=True)

exp.diag_table = diag

# define namelist values as python dictionary
# wrapped as a namelist object.
namelist = Namelist({
    'main_nml': {
        'dt_atmos': 300.,#128./1.,
        'days': 100./1.,
        'calendar': 'no_calendar',
        'current_date': [1,1,1,0,0,0]
    },

    'atmosphere_nml': {
        'idealized_moist_model': False  # False for Newtonian Cooling.  True for Isca/Frierson
    },


    'astronomy_nml':{
        'obliq':0.0, 
        'ecc':0.0,
    }, 

    'constants_nml':{
        'omega'   : omega,
        'orbital_period':2*3.14159265358979323846/omega,
        #'grav'    : grav, 
        'radius'  : radius, 
        #'rdgas'   : Rd, 
        #'PSTD'    : ps * 10., 
        #'PSTD_MKS': ps, 
        #'kappa'   : Rd / cp,
    },
        

    # 'spectral_dynamics_nml': {
    #     'damping_order'           : 4,                      # default: 2
    #     'water_correction_limit'  : 200.e2,                 # default: 0
    #     'reference_sea_level_press': ps,                  # default: 101325
    #     'valid_range_t'           : [100., 800.],           # default: (100, 500)
    #     'initial_sphum'           : 0.0,                  # default: 0
    #     'vert_coord_option'       : 'even_sigma',         # default: 'even_sigma'
    #     'scale_heights': 6.0,
    #     'exponent': 7.5,
    #     'surf_res': 0.5
    # },

    # 'spectral_dynamics_nml':{
    #     'damping_order':4, 
    #     #'damping_coeff':1/86400., 
    #     #'damping_option':'exponential_cutoff', 
    #     #'cutoff_wn':15, 
    #     'do_water_correction':False, 
    #     'water_correction_limit'
    #     'reference_sea_level_press':ps, 
    #     'initial_sphum':0.0, 
    #     'valid_range_t':[100., 800.], 
    #     'robert_coeff':0.03, 
    #     'num_levels':20, 
    #     'vert_coord_option':'even_sigma', 
    #     #'num_levels':25, 
    #     #'vert_coord_option':'hybrid', 
    #     #'surf_res':0.05, 
    #     #'scale_heights':5.0, 
    #     #'exponent':3.0, 
    #     'lat_max': 64,#*4,#*2, 
    #     'lon_max': 128,#*4,#*2, 
    #     'num_fourier':42,#170,#85, 
    #     'num_spherical':42+1,#170+1,#85+1  
    #     #'raw_filter_coeff':0.53
    # },

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
        'vert_coord_option'       : 'uneven_sigma',         # default: 'even_sigma'
        'scale_heights': 6.0,
        'exponent': 7.5,
        'surf_res': 0.5,
        'num_levels':25, 
        'lat_max': 64,
        'lon_max': 128,
        'num_fourier':42,
        'num_spherical':42+1
    },

    # configure the relaxation profile
    'hs_forcing_nml': {
        'equilibrium_t_option':'EXOPLANET_TL',
        't_zero': 315.,    # temperature at reference pressure at equator (default 315K)
        't_strat': 200.,   # stratosphere temperature (default 200K)
        'delh': 115.,       # equator-pole temp gradient (default 60K)
        'delv': 10.,       # lapse rate (default 10K)
        'eps': 0.,         # stratospheric latitudinal variation (default 0K)
        'sigma_b': 0.7,    # boundary layer friction height (default p/ps = sigma = 0.7)

        # negative sign is a flag indicating that the units are days
        'ka':   -20./1.,#*32*np.sqrt(2),      # Constant Newtonian cooling timescale (default 40 days)
        'ks':    -5./1.,#*32*np.sqrt(2),      # Boundary layer dependent cooling timescale (default 4 days)
        'kf':   -1./1.,       # BL momentum frictional timescale (default 1 days)

        'do_conserve_energy':   True,  # convert dissipated momentum into heat (default True)
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
#exp.set_resolution(*RESOLUTION)

#oms = 2**np.linspace(1.5,6,10)
om = 4. 
rads = np.array([1./2.,1./2.,1./4.])
names = ['4_1tau','2_1tau','4_1tau']#,'1_2tau','1_4tau']
radnames = ['1_2rad', '1_2rad', '1_4rad']
taus  = np.array([4.,2.,4.])#,0.5,0.25])
steps = np.array([64.,64.,64.])#,64.,64.])#32.,32.,32.,64.,128.,256.,512.,512.])

#Lets do a run!
if __name__ == '__main__':

    for n, tau in enumerate(taus):
        omega_exp = omega_e / om
        rad_exp = radius * rads[n]
        ka_exp = -20. * tau 
        ks_exp = -5. * tau 
        step_exp = steps[n]
        run_exp = exp.derive('tl_'+radnames[n]+'_'+'1_4om_'+names[n]+'_rad_diff10')
        run_exp.namelist['main_nml']['dt_atmos'] = step_exp
        run_exp.namelist['constants_nml']['omega'] = omega_exp
        run_exp.namelist['constants_nml']['radius'] = rad_exp 
        run_exp.namelist['hs_forcing_nml']['ka'] = ka_exp 
        run_exp.namelist['hs_forcing_nml']['ks'] = ks_exp
        run_exp.run(1, num_cores=NCORES,use_restart=False,overwrite_data=False)
        for i in range(2, 21):
            run_exp.run(i, num_cores=NCORES, overwrite_data=False)  # use the restart i-1 by default
