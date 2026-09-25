import os

import numpy as np

from isca import ShallowCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES = 16
base_dir = os.path.dirname(os.path.realpath(__file__))
# a CodeBase can be a directory on the computer,
# useful for iterative development
cb = ShallowCodeBase.from_directory(GFDL_BASE)

# or it can point to a specific git repo and commit id.
# This method should ensure future, independent, reproducibility of results.
# cb = DryCodeBase.from_repo(repo='https://github.com/isca/isca', commit='isca1.1')

# compilation depends on computer specific settings.  The $GFDL_ENV
# environment variable is used to determine which `$GFDL_BASE/src/extra/env` file
# is used to load the correct compilers.  The env file is always loaded from
# $GFDL_BASE and not the checked out git repo.

# create an Experiment object to handle the configuration of model parameters
# and output diagnostics
exp = Experiment('shallow_pb_showman_like', codebase=cb)

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_out', 1, 'days', time_units='days')
diag.add_field('shallow_diagnostics', 'ucomp', time_avg=False)
diag.add_field('shallow_diagnostics', 'vcomp', time_avg=False)
diag.add_field('shallow_diagnostics', 'vor', time_avg=False)
diag.add_field('shallow_diagnostics', 'div', time_avg=False)
diag.add_field('shallow_diagnostics', 'h', time_avg=False)
diag.add_field('shallow_diagnostics', 'pv_corrected', time_avg=False)
diag.add_field('shallow_diagnostics', 'stream', time_avg=False)
diag.add_field('shallow_diagnostics', 'du_dt_mass', time_avg=False)
diag.add_field('shallow_diagnostics', 'dv_dt_mass', time_avg=False)
diag.add_field('shallow_diagnostics', 'h_eq', time_avg=False)
diag.add_field('shallow_diagnostics', 'du_dt_drag', time_avg=False)
diag.add_field('shallow_diagnostics', 'dv_dt_drag', time_avg=False)

exp.diag_table = diag

#Empty the run directory ready to run
exp.clear_rundir()

#Define values for the 'core' namelist
exp.namelist = namelist = Namelist({
  'main_nml':{
   'days'   : 20,
   'hours'  : 0,
   'minutes': 0,
   'seconds': 0,
   'dt_atmos': 600,
   'calendar': 'no_calendar',
  },

 'atmosphere_nml':{
   'print_interval': 86400,
    },

'fms_io_nml':{
   'threading_write' :'single',
   'fileset_write': 'single'
    },

 'fms_nml':{
   'print_memory_usage':True,
   'domains_stack_size': 200000,
    },

 'shallow_dynamics_nml':{
   'num_lon'             : 512, # 256
   'num_lat'             : 256, # 128
   'num_fourier'         : 170, # 85
   'num_spherical'       : 171, # 86
   'fourier_inc'         : 1,
   'damping_option'      : 'resolution_dependent',
   'damping_order'       : 4,
   'damping_coeff'       : (1./0.1/86400.), 
   'h_0'                 : 4.e6,
   'grid_tracer'         : False,
   'spec_tracer'         : False,
   'robert_coeff'        : 0.04,
   'robert_coeff_tracer' : 0.04,
    },

 'shallow_physics_nml': {
   'h_eq_option'      : 'perez_becker',  # default 'legacy' 
   'do_zero_mean_h_eq': True,
   'h_lon'            : 180.0,                # default 90.0 (!) — substellar longitude
   'del_h'            : 0.1,                # default 0.0 — must be < 1
   'h_0'              : 4.e6,               # must equal shallow_dynamics_nml h_0
   'therm_damp_time'  : 0.1*86400.,         # negative = days, positive = seconds, 0 disables
   'fric_damp_time'   : 10.*86400.,         # same convention; 0.0 for no drag
   'do_mass_exchange' : True,               # default False
   },
  
  'constants_nml': { 
    'radius': 8.2e7, 
    'omega': 3.2e-5, 
  }, 

})

#Lets do a run!
if __name__=="__main__":
    cb.compile()  # compile the source code to working directory $GFDL_WORK/codebase

    exp.run(1, use_restart=False, num_cores=NCORES)
    for i in range(2,6):
        exp.run(i, num_cores=NCORES)