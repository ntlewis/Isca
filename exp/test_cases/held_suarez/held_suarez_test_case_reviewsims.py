import numpy as np

from isca import GreyCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES = 32
RESOLUTION = 'T85', 15  # T42 horizontal resolution, 25 levels in pressure

# a CodeBase can be a directory on the computer,
# useful for iterative development
cb = GreyCodeBase.from_directory(GFDL_BASE)

# or it can point to a specific git repo and commit id.
# This method should ensure future, independent, reproducibility of results.
# cb = DryCodeBase.from_repo(repo='https://github.com/isca/isca', commit='isca1.1')

# compilation depends on computer specific settings.  The $GFDL_ENV
# environment variable is used to determine which `$GFDL_BASE/src/extra/env` file
# is used to load the correct compilers.  The env file is always loaded from
# $GFDL_BASE and not the checked out git repo.

#cb.compile()  # compile the source code to working directory $GFDL_WORK/codebase

# create an Experiment object to handle the configuration of model parameters
# and output diagnostics

exp_name = 'held_suarez_review'
exp = Experiment(exp_name, codebase=cb)

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 100, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'height', time_avg=True)
diag.add_field('dynamics', 'omega', time_avg=True)

exp.diag_table = diag

# define namelist values as python dictionary
# wrapped as a namelist object.
namelist = Namelist({
    'main_nml': {
        'dt_atmos': 120,
        'days': 100,
        'calendar': 'no_calendar',
        'current_date': [1,1,1,0,0,0]
    },

    'atmosphere_nml': {
        'idealized_moist_model': False  # False for Newtonian Cooling.  True for Isca/Frierson
    },

    'spectral_dynamics_nml': {
        'damping_order'           : 4,                      # default: 2
        'damping_coeff'           : 10./86400., 
        'water_correction_limit'  : 200.e2,                 # default: 0
        'reference_sea_level_press': 1.0e5,                  # default: 101325
        'valid_range_t'           : [100., 800.],           # default: (100, 500)
        'initial_sphum'           : 0.0,                  # default: 0
        'vert_coord_option'       : 'even_sigma',         # default: 'even_sigma'
        'scale_heights': 6.0,
        'exponent': 7.5,
        'surf_res': 0.5
    },

    # configure the relaxation profile
    'hs_forcing_nml': {
        't_zero': 315.,    # temperature at reference pressure at equator (default 315K)
        't_strat': 200.,   # stratosphere temperature (default 200K)
        'delh': 60.,       # equator-pole temp gradient (default 60K)
        'delv': 10.,       # lapse rate (default 10K)
        'eps': 0.,         # stratospheric latitudinal variation (default 0K)
        'sigma_b': 0.7,    # boundary layer friction height (default p/ps = sigma = 0.7)

        # negative sign is a flag indicating that the units are days
        'ka':   -40.,      # Constant Newtonian cooling timescale (default 40 days)
        'ks':    -4.,      # Boundary layer dependent cooling timescale (default 4 days)
        'kf':   -1.,       # BL momentum frictional timescale (default 1 days)

        'do_conserve_energy':   True,  # convert dissipated momentum into heat (default True)
    },

    'constants_nml':{
        'radius':6.371e6
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
exp.set_resolution(*RESOLUTION)

rads = np.array([0.125])
rad_names = ['1_8rad']
taurs = np.array([1., 0.1])
taur_names = ['1_1tau_rad', '1_10tau_rad']
tauds = np.array([1., 0.1])
taud_names = ['1_1tau_d', '1_10tau_d']

#Lets do a run!
if __name__ == '__main__':

    for ridx, rad in enumerate(rads):
        for tridx, taur in enumerate(taurs):
            for tdidx, taud in enumerate(tauds):
                run_name = 'hsreview_'+rad_names[ridx]+'_1_1om_'+taur_names[tridx]+'_'+taud_names[tdidx]
                run_exp = exp.derive(run_name)
                if rad == 2:
                    RESOLUTION = 'T85', 15
                    run_exp.set_resolution(*RESOLUTION)
                else:
                    RESOLUTION = 'T42', 15
                    run_exp.set_resolution(*RESOLUTION)
                run_exp.namelist['constants_nml']['radius'] = 6.371e6 * rad 
                run_exp.namelist['hs_forcing_nml']['ka'] = -40. * taur
                run_exp.namelist['hs_forcing_nml']['ks'] = -4.  * taur
                run_exp.namelist['hs_forcing_nml']['kf'] = -1. * taud 
                run_exp.namelist['spectral_dynamics_nml']['damping_coeff'] = 10./86400. / taud 
                run_exp.run(1, num_cores=NCORES, use_restart=False, overwrite_data=True)
                for i in range(2,11): 
                    if i == 10:
                        diag = DiagTable()
                        diag.add_file('atmos_monthly', 1, 'days', time_units='days')
                        #Tell model which diagnostics to write
                        diag.add_field('dynamics', 'ps', time_avg=True)
                        diag.add_field('dynamics', 'bk')
                        diag.add_field('dynamics', 'pk')
                        diag.add_field('dynamics', 'ucomp', time_avg=True)
                        diag.add_field('dynamics', 'vcomp', time_avg=True)
                        diag.add_field('dynamics', 'temp', time_avg=True)
                        diag.add_field('dynamics', 'height', time_avg=True)
                        diag.add_field('dynamics', 'omega', time_avg=True)
                        run_exp.diag_table = diag
                    run_exp.run(i, num_cores=NCORES, overwrite_data=True)
                    
