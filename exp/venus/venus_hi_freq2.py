import numpy as np

from isca import DryCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES = 16
RESOLUTION = 'T21', 31  # change here from T21 horizontal resolution, 31 levels in pressure

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

exp_name = 'venus_hi_freq3'
exp = Experiment(exp_name, codebase=cb)

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_out', 6, 'hours', time_units='days')

#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=False)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('dynamics', 'ucomp', time_avg=False)
diag.add_field('dynamics', 'vcomp', time_avg=False)
diag.add_field('dynamics', 'temp', time_avg=False)
diag.add_field('dynamics', 'vor', time_avg=False)
diag.add_field('dynamics', 'div', time_avg=False)
diag.add_field('hs_forcing', 'teq', time_avg=False)
diag.add_field('dynamics', 'height', time_avg=False)
diag.add_field('dynamics', 'height_half', time_avg=False)
diag.add_field('dynamics', 'pres_full', time_avg=False)
diag.add_field('dynamics', 'pres_half', time_avg=False)

# diagnostics for angular momentum budget 
diag.add_field('dynamics', 'dt_ug_damp', time_avg=False) 
diag.add_field('hs_forcing', 'udt_rdamp', time_avg=False)
diag.add_field('dynamics', 'dt_ug_dyn', time_avg=False)
diag.add_field('dynamics', 'dt_psg_tot', time_avg=False)
diag.add_field('dynamics', 'dt_ug', time_avg=False)
# diagnostics for EP flux 
diag.add_field('dynamics', 'up_vp', time_avg=False)
diag.add_field('dynamics', 'theta', time_avg=False)
diag.add_field('dynamics', 'up_omegap', time_avg=False)
diag.add_field('dynamics', 'vp_thetap', time_avg=False)



exp.diag_table = diag

# define namelist values as python dictionary
# wrapped as a namelist object.

etas = np.array([2.265e-6, 2.925e-6, 4.500e-6, 7.950e-6, 1.525e-5, 
3.070e-5, 6.355e-5, 1.324e-4, 2.710e-4, 5.400e-4, 1.043e-3, 1.934e-3, 3.454e-3, 
5.929e-3, 9.814e-3, 1.569e-2, 2.430e-2, 3.657e-2, 5.361e-2, 7.672e-2, 1.074e-1, 
1.472e-1, 1.979e-1, 2.610e-1, 3.373e-1, 4.273e-1, 5.299e-1, 6.420e-1, 7.577e-1, 
8.679e-1, 9.602e-1])
#etas = np.array([0.00022, 0.00029, 0.00045, 0.00078, 0.0015, 0.0030, 0.0062, 0.0130, 0.0266, 0.0532, 0.1029, 0.1912, 0.3419, 0.5878, 0.9739, 1.5584, 2.4162, 3.6385, 5.3373, 7.6420, 10.703, 14.678, 19.738, 26.035, 33.660, 42.658, 52.918, 64.136, 75.715, 86.746, 95.994]) * 1.e-2
etas = (etas[:-1] + etas[1:])/2.

half_etas = np.concatenate(([0.],etas,[1.])) 
half_etas = half_etas.tolist()

namelist = Namelist({
    'main_nml': {
        'dt_atmos': 150, # change here from 150
        'days': 10,
        'current_date' : [1,1,1,0,0,0],
        'calendar': 'no_calendar',
    },

    'atmosphere_nml': {
        'idealized_moist_model': False  # Use HS Forcing, no WV, no surface. 
    },

    'spectral_dynamics_nml': {                                                      # NTL: (**) indicates action required 
        'damping_order'           : 4,                      # default: 2            # NTL: As suggested in Lee and Richardson (2010)
        'damping_coeff'           : 4e-6,                   # default: 1.157...e-4  # NTL: Diffusion used by LR 2010. 
        'water_correction_limit'  : 200.e2,                 # default: 0            # NTL: What does this do?                          **
        'reference_sea_level_press': 9.2e6,                 # default: 101325       # NTL: Change to Venus surface pressure 
        'valid_range_t'           : [50., 850.],            # default: (100, 500)   # NTL: extend range for Venus        
        'initial_sphum'           : 0.0,                    # default: 0            # NTL: Start with no water
        'vert_coord_option'       : 'input',                # default: 'even_sigma' # NTL: change this to input    
    },

    # CONVERT LEE INPUT PRESSURE INTO INPUT SIGMA... 
    'vert_coordinate_nml': {
       'bk': half_etas,
       'pk': np.zeros(len(half_etas)).tolist(),
      },

    # configure the relaxation profile -- for the time being, read this from input, use what Chris and Joao used. 
    'hs_forcing_nml': {
        'equilibrium_t_option' : 'venus',
        'venus_model'          : True, 
        'do_conserve_energy':   True,  # convert dissipated momentum into heat (default True)
    },

    'constants_nml': {
        'omega': 2.99e-7, 
        'grav': 8.87,
        'radius': 6.0518e6,
        'rdgas': 188.,
        'PSTD':92.0e6,
        'PSTD_MKS':92.0e5,
        'kappa': 188. / 850.1,
        'wtmair':43.45
        #'CP_AIR': 900.,

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

#Lets do a run!
if __name__ == '__main__':
    exp.run(1, num_cores=NCORES, restart_file='/network/group/aopp/planetary/PLR016_LEWIS_VENUS/isca_data/venus_default5/restarts/res0515.tar.gz')
    for i in range(2, 201): # change here from 601
        exp.run(i, num_cores=NCORES)  # use the restart i-1 by default
    #exp.run(487, num_cores=NCORES, use_restart=True, overwrite_data=True)
    #for i in range(488, 501):
        #exp.run(i, num_cores=NCORES, overwrite_data=True)
