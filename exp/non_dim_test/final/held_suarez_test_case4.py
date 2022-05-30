import numpy as np

from isca import DryCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES = 48
RESOLUTION = 'T63', 35  # T42 horizontal resolution, 25 levels in pressure

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

exp_name = 'held_suarez_default'
exp = Experiment(exp_name, codebase=cb)

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_out', 30, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'vor', time_avg=True)
diag.add_field('dynamics', 'div', time_avg=True)

exp.diag_table = diag

# define namelist values as python dictionary
# wrapped as a namelist object.
namelist = Namelist({
    'main_nml': {
        'dt_atmos': 150,
        'days': 360,
        'calendar': 'thirty_day',
        'current_date': [2000,1,1,0,0,0]
    },

    'atmosphere_nml': {
        'idealized_moist_model': False  # False for Newtonian Cooling.  True for Isca/Frierson
    },

    'spectral_dynamics_nml': {
        'damping_order'           : 4,                      # default: 2
        'water_correction_limit'  : 200.e2,                 # default: 0
        'reference_sea_level_press': 1.0e5,                  # default: 101325
        'valid_range_t'           : [100., 800.],           # default: (100, 500)
        'initial_sphum'           : 0.0,                  # default: 0
        'vert_coord_option'       : 'uneven_sigma',         # default: 'even_sigma'
        'scale_heights': 6.0,
        'exponent': 7.5,
        'surf_res': 0.5
    },

    # configure the relaxation profile
    'hs_forcing_nml': {
        'equilibrium_t_option':'Wang_Read',
        't_zero': 288.,    # temperature at reference pressure at equator (default 315K)
        #'t_strat': 200.,   # stratosphere temperature (default 200K)
        'delh': 60.,       # equator-pole temp gradient (default 60K)
        #'delv': 10.,       # lapse rate (default 10K)
        'smoothing_T': 2., 
        'z_trop': 12000.,
        #'eps': 0.,         # stratospheric latitudinal variation (default 0K)
        'sigma_b': 0.8,    # boundary layer friction height (default p/ps = sigma = 0.7)

        # negative sign is a flag indicating that the units are days
        'ka':   -30.,      # Constant Newtonian cooling timescale (default 40 days)
        'ks':    -2.5,      # Boundary layer dependent cooling timescale (default 4 days)
        'kf':   -0.6,       # BL momentum frictional timescale (default 1 days)

        'do_conserve_energy':   True,  # convert dissipated momentum into heat (default True)
        'sponge':True, 
        'sp_coeff': -0.5, 
        'sigma_sp': 10. / 1000.
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
    }, 
    'constants_nml':{
        'radius': 6400.e3, 
        'grav':9.8, 
        'omega':7.27e-5
    }
})

exp.namelist = namelist
exp.set_resolution(*RESOLUTION)

#Lets do a run!
if __name__ == '__main__':
    rads = [23.]#, 8., 23.]
    earth_rad = 6400.e3
    earth_omega = 7.27e-5
    n_exp = 1
    exps_dict={}
    exps_dict[0] = {'rad_factor':16, 'omg_factor':16,'dt':150}
    for n in range(n_exp):
        exp_dict = exps_dict[n]
        rad_factor = exp_dict['rad_factor']
        omg_factor = exp_dict['omg_factor']
        dt         = exp_dict['dt']
        nyears = int(7 * omg_factor)
        run_exp = exp.derive('1_wang_a_'+str(rad_factor)+'_omg_'+str(omg_factor))
        run_exp.namelist['constants_nml']['radius'] = earth_rad / rad_factor
        run_exp.namelist['constants_nml']['omega'] = earth_omega / omg_factor
        run_exp.namelist['spectral_dynamics_nml']['damping_coeff'] = 1 / (0.25 * 86400.) / omg_factor
        run_exp.namelist['hs_forcing_nml']['ka'] = omg_factor * -30.
        run_exp.namelist['hs_forcing_nml']['ks'] = omg_factor * -2.5
        run_exp.namelist['hs_forcing_nml']['kf'] = omg_factor * -0.6
        run_exp.namelist['hs_forcing_nml']['sp_coeff'] = omg_factor * -0.5
        run_exp.namelist['main_nml']['dt_atmos'] = dt
        run_exp.run(1, num_cores=NCORES, use_restart=False, overwrite_data=False)
        for i in range(2, nyears+1):
            if i >= nyears - 2.5:
                #Tell model how to write diagnostics
                diag = DiagTable()
                diag.add_file('atmos_out', 6, 'hours', time_units='days')
                #Tell model which diagnostics to write
                diag.add_field('dynamics', 'ps', time_avg=True)
                diag.add_field('dynamics', 'bk')
                diag.add_field('dynamics', 'pk')
                diag.add_field('dynamics', 'ucomp', time_avg=True)
                diag.add_field('dynamics', 'vcomp', time_avg=True)
                diag.add_field('dynamics', 'omega', time_avg=True)
                diag.add_field('dynamics', 'temp', time_avg=True)
                diag.add_field('dynamics', 'vor', time_avg=True)
                diag.add_field('dynamics', 'div', time_avg=True)
                diag.add_field('hs_forcing', 'teq', time_avg=True)
                diag.add_field('dynamics', 'height', time_avg=True)
                diag.add_field('dynamics', 'height_half', time_avg=True)
                diag.add_field('dynamics', 'pres_full', time_avg=True)
                diag.add_field('dynamics', 'pres_half', time_avg=True)
                # diagnostics for angular momentum budget 
                diag.add_field('dynamics', 'dt_ug_damp', time_avg=True) 
                diag.add_field('hs_forcing', 'udt_rdamp', time_avg=True)
                diag.add_field('dynamics', 'dt_ug_dyn', time_avg=True)
                diag.add_field('dynamics', 'dt_psg_tot', time_avg=True)
                diag.add_field('dynamics', 'dt_ug', time_avg=True)
                # diagnostics for EP flux 
                diag.add_field('dynamics', 'up_vp', time_avg=True)
                diag.add_field('dynamics', 'theta', time_avg=True)
                diag.add_field('dynamics', 'up_omegap', time_avg=True)
                diag.add_field('dynamics', 'vp_thetap', time_avg=True)
                # add diag_table to experiment 
                run_exp.diag_table = diag
            run_exp.run(i, num_cores=NCORES, overwrite_data=False)  # use the restart i-1 by default
