"""
This python script configures builds the idealised Venus GCM (spectral model) 
described in Lee and Richardson (2010) using Isca (Vallis et al., 2018).  

Radiative transfer is parameterised as Newtonian cooling that relaxes to a 
temperature field derived from VIRA (the Venus international reference atmosphere),
following Lee et al. (2007). 

Linear rayleigh drag is applied at the lowest model level to parameterise the 
boundary layer, and to the eddies at the model top (top three levels) as a sponge
layer. This approach follows Lee et al. (2007)

The boundary layer and sponge parametrisations are called by setting venus_model = True 
in hs_forcing_nml. The temperature relaxation field is configured here to be passed to 
the model as an input file, 'venus_relaxation_profile.nc', but it can also be configured 
by setting equilibrium_t_option = 'venus'. 

Finally, orbital, atmospheric, and planetary parameters are set to be those of Venus following 
Lee (2006). 

References:

1) Lee, C. (2006) Modelling of the Atmosphere of Venus, D.Phil. Thesis, University of Oxford, 
    Oxford, UK

2) Lee, C., S. R. Lewis, and P. L. Read (2007), Superrotation in a Venus GCM, J. Geophys. Res., 
    112, E04S11

3) Lee, C., and M. I. Richardson (2010), A general circulation model ensemble study of the atmospheric 
    circulation of Venus, J. Geophys. Res., 115, E04002 

4) Vallis, G. K., and Coauthors (2018), Isca, v1.0: A framework for the global modelling of the atmospheres 
    of Earth and other planets at varying levels of complexity. Geophys. Model Dev., 11, 843–859

Please direct questions to Neil Lewis at neil.lewis@physics.ox.ac.uk 
NTL 13/11/2018
"""

# import required python modules 
import numpy as np
import os
from isca import DryCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE


# locate and compile source code, configure expeiment object. 
cb = DryCodeBase.from_directory(GFDL_BASE) # tell model where code is
cb.compile()                               # compile the source code to working directory $GFDL_WORK/codebase
exp_name = 'venus_test_case'               # name experiment 
exp = Experiment(exp_name, codebase=cb)    # create an Experiment object to handle the configuration of model parameters
                                           # and output diagnostics

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_out', 75, 'days', time_units='days')
#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
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
exp.diag_table = diag



# load temperature relaxation file, and construct model pressure levels.
exp.inputfiles = [os.path.join(GFDL_BASE,'input/venus_input_files/venus_relaxation_profile.nc')]
# define vertical levels as used by Lee and Richardson (2010, hereafter LR2010)
etas = np.array([2.265e-6, 2.925e-6, 4.500e-6, 7.950e-6, 1.525e-5, 
3.070e-5, 6.355e-5, 1.324e-4, 2.710e-4, 5.400e-4, 1.043e-3, 1.934e-3, 3.454e-3, 
5.929e-3, 9.814e-3, 1.569e-2, 2.430e-2, 3.657e-2, 5.361e-2, 7.672e-2, 1.074e-1, 
1.472e-1, 1.979e-1, 2.610e-1, 3.373e-1, 4.273e-1, 5.299e-1, 6.420e-1, 7.577e-1, 
8.679e-1, 9.602e-1])
central_half_etas = (etas[:-1] + etas[1:])/2.
half_etas = np.concatenate(([0.],central_half_etas,[1.])) 
half_etas = half_etas.tolist()

### configure experiment resolution and number of cores to use 
NCORES = 16
RESOLUTION = 'T21', 31

# define namelist values 
namelist = Namelist({
    'main_nml': {
        'dt_atmos': 150, 
        'days': 75,
        'current_date' : [1,1,1,0,0,0],
        'calendar': 'no_calendar',
    },

    'atmosphere_nml': {
        'idealized_moist_model': False  # Use HS Forcing, no WV, no surface. 
    },

    'spectral_dynamics_nml': {                                                      
        'damping_order'           : 4,                      # As suggested in LR2010,             default: 2     
        'damping_coeff'           : 4e-6,                   # Diffusion used by LR2010,           default: 1.157...e-4  
        'water_correction_limit'  : 200.e2,                 #                                     default: 0            
        'reference_sea_level_press': 9.2e6,                 # Change to Venus surface pressure,   default: 101325         
        'valid_range_t'           : [50., 850.],            # extend valid temp. range for Venus, default: (100, 500)          
        'initial_sphum'           : 0.0,                    # Start with no water,                default: 0           
        'vert_coord_option'       : 'input',                # Use pressure levels from L2007,     default: 'even_sigma' 
    },

    
    'vert_coordinate_nml': {
       'bk': half_etas, # Pressure levels constructed above, as defined by LR2007
       'pk': np.zeros_like(half_etas).tolist(),
      },

    # configure the relaxation profile 
    'hs_forcing_nml': {
        'equilibrium_t_option' : 'from_file',
        'equilibrium_t_file'   : 'venus_relaxation_profile', # Relaxation profile as defined by L2007       
        'venus_model'          : True,                       # Use rayleigh boundary layer and sponge (eddies only)
        'do_conserve_energy'   : True,                       # convert dissipated momentum into heat (default True)
    },

    'constants_nml': { # set constants to be those of Venus 
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

# Run the model. NB: This is set to run for a LONG time, i.e. until the model has spun-up. 
if __name__ == '__main__':
    exp.run(1, num_cores=NCORES, use_restart=False)
    for i in range(2, 601): # run for a long time, required to completely spin-up model. 
        exp.run(i, num_cores=NCORES)  
