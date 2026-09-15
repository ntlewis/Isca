"""
Barotropic vorticity equation experiment configured for forced-dissipative
turbulence on the sphere, following the approach of:

  Scott, R. K. and Polvani, L. M., 2007: Forced-dissipative shallow-water
  turbulence on the sphere and the atmospheric circulation of the giant planets.
  J. Atmos. Sci., 64, 3158-3176.

This script uses the barotropic (non-divergent) vorticity equation rather than
the full shallow-water system, which corresponds to the infinite-deformation-
radius limit of S&P. The key feature is white-in-time (delta-correlated)
isotropic stochastic forcing in a narrow band of total wavenumbers, as opposed
to the Ornstein-Uhlenbeck (Markovian) forcing used in the original Vallis et al.
(2004) experiment.

FORCING CONFIGURATION (stirring_nml):
  do_white_forcing  = True   -> white-in-time forcing (S&P), delta-correlated in time
  do_localize       = False  -> isotropic: no latitude/longitude mask
  zonal_forcing_min = 0      -> force all zonal wavenumbers m >= 1
  n_total_forcing_min/max    -> sets the forced total-wavenumber band
                                forced wavenumbers L satisfy n_min < L < n_max

AMPLITUDE / ENERGY INJECTION RATE:
  Specify epsilon_0 (m^2/s^3) in stirring_nml directly. The Fortran code
  computes the stirring amplitude internally from:

    amplitude = sqrt( 2 * epsilon_0 / (dt * a^2 * S) )

  where:
    S   = sum_{L = n_total_forcing_min+1}^{n_total_forcing_max-1} 1/(L+1)
    dt  = dt_atmos (model timestep, s)
    a   = planetary radius (m)  [from constants_mod]

  This formula is calibrated so that the kinetic energy injection rate
  equals epsilon_0 (m^2/s^3) in the GFDL leapfrog scheme, where the
  leapfrog step is 2*dt_atmos while stirring_init receives dt_atmos.

  To convert from S&P's non-dimensional epsilon* (non-dimensionalised by
  a^2 * Omega^3, where Omega is the planetary rotation rate in rad/s):
    epsilon_0 = epsilon* * a^2 * Omega^3

DISSIPATION (barotropic_dynamics_nml):
  damping_coeff_r  -> linear (Rayleigh) drag rate (s^-1), acts on all wavenumbers.
                      Corresponds to alpha in S&P. A timescale of 10-50 days is typical.
  damping_order    -> order p of hyperdiffusion operator (-1)^p * nabla^(2p).
                      p=4 gives nabla^8 (recommended for scale-separated turbulence).
  damping_coeff    -> hyperdiffusion rate for the highest retained mode.
                      Set via damping_option='resolution_dependent'.

TYPICAL PARAMETER CHOICES:
  To reproduce S&P-style jets, vary:
    - The forcing wavenumber band (larger L_f -> smaller injection scale)
    - The drag timescale 1/damping_coeff_r (larger -> more energy, stronger jets)
    - epsilon_0 (controls energy injection rate and equilibrium energy level)
  The Rhines scale L_beta ~ sqrt(U/beta) sets the expected jet width.
"""

import os
import numpy as np
from isca import BarotropicCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES = 32
base_dir = os.path.dirname(os.path.realpath(__file__))
cb = BarotropicCodeBase.from_directory(GFDL_BASE)

exp = Experiment('barotropic_stochastic_white_forcing', codebase=cb)

# -------------------------------------------------------------------------
# Diagnostics
# -------------------------------------------------------------------------
diag = DiagTable()
diag.add_file('atmos_monthly', 1, 'days', time_units='days')

diag.add_field('barotropic_diagnostics', 'ucomp',    time_avg=False)
diag.add_field('barotropic_diagnostics', 'vcomp',    time_avg=False)
diag.add_field('barotropic_diagnostics', 'vor',      time_avg=False)
#diag.add_field('barotropic_diagnostics', 'pv',       time_avg=True)
diag.add_field('barotropic_diagnostics', 'stream',   time_avg=False)
#diag.add_field('barotropic_diagnostics', 'trs',      time_avg=True)
#diag.add_field('barotropic_diagnostics', 'tr',       time_avg=True)
#diag.add_field('barotropic_diagnostics', 'eddy_vor', time_avg=True)
#diag.add_field('barotropic_diagnostics', 'delta_u',  time_avg=True)
#diag.add_field('stirring_mod',           'stirring',     time_avg=True)
#diag.add_field('stirring_mod',           'stirring_amp', time_avg=True)
#diag.add_field('stirring_mod',           'stirring_sqr', time_avg=True)

exp.diag_table = diag
exp.clear_rundir()

# -------------------------------------------------------------------------
# Namelist
# -------------------------------------------------------------------------

# --- FORCING BAND ---
# Forced total wavenumbers l satisfy: n_total_forcing_min < l < n_total_forcing_max
# i.e., l in [n_total_forcing_min+1, n_total_forcing_max-1].
# Central wavenumber l_f = (n_total_forcing_min + n_total_forcing_max) / 2.
# Example below: l_f = 20, band l = 18..22 (width = 5 wavenumbers).
# For a higher forcing wavenumber (smaller injection scale), increase both values.
N_FORCING_MIN = 45   # l_f = (17 + 23) / 2 = 20,  band l = 18..22
N_FORCING_MAX = 51

# --- ENERGY INJECTION RATE ---
# epsilon_0 in m^2/s^3. Passed directly to stirring_nml.
# The Fortran code computes amplitude internally.
# To convert from S&P non-dimensional epsilon* (non-dimensionalised by a^2 * Omega^3):
#   epsilon_0 = epsilon* * a^2 * Omega^3
#   e.g. epsilon* = 1e-4, a = 6.371e6 m, Omega = 7.292e-5 rad/s:
#   epsilon_0 = 1e-4 * (6.371e6)^2 * (7.292e-5)^3 = ...
epsilon_0 = 5 * 4.e-9    # m^2/s^3



# --- LINEAR DRAG ---
# damping_coeff_r is the Rayleigh drag rate (s^-1).
# 1/(20 days) = 5.79e-7 s^-1.  Longer timescale -> stronger jets.
DRAG_TIMESCALE_DAYS = 1.e4
damping_coeff_r = 5.0 / (DRAG_TIMESCALE_DAYS * 86400.0)  # s^-1

exp.namelist = namelist = Namelist({

    'main_nml': {
        'days'    : 1000,#1000,
        'hours'   : 0,
        'minutes' : 0,
        'seconds' : 0,
        'dt_atmos': 1200,
        'calendar': 'no_calendar',
        'print_memuse':False,
    },

    'atmosphere_nml': {
        'print_interval': 86400,
    },

    'fms_io_nml': {
        'threading_write': 'single',
        'fileset_write'  : 'single',
    },

    'fms_nml': {
        'print_memory_usage': True,
        'domains_stack_size': 200000,
    },

    'barotropic_dynamics_nml': {
        'triang_trunc'   : True,
        'num_lat'        : 288,
        'num_lon'        : 576,
        'num_fourier'    : 191,
        'num_spherical'  : 192,
        'fourier_inc'    : 1,

        # Hyperdiffusion: damping_order=4 gives nabla^8, appropriate for
        # turbulence experiments where clean scale separation is required.
        # Use damping_order=2 (nabla^4) if preferred.
        'damping_option' : 'resolution_dependent',
        'damping_order'  : 4,          # nabla^8 hyperdiffusion
        'damping_coeff'  : 1. / (86400. / 10.),    # rate for highest retained mode (s^-1)

        # Linear Rayleigh drag: acts uniformly on all wavenumbers.
        # Corresponds to alpha in S&P. Tune to control jet strength.
        'damping_coeff_r': damping_coeff_r,
        'grid_tracer'  : False,   # disable tracers 
        'spec_tracer'  : False,

        # Start from rest (no initial jets or eddies).
        # The turbulence will spin up from the stochastic forcing alone.
        'initial_zonal_wind': 'zero',
        'zeta_0'            : 0.0,    # no initial eddy perturbation
        'm_0'               : 4,      # (unused when zeta_0=0)
        'eddy_lat'          : 45.0,
        'eddy_width'        : 10.0,

        'robert_coeff'  : 0.04,
    },

    'barotropic_physics_nml': {},

    'stirring_nml': {
        # White-in-time S&P forcing
        'do_white_forcing'   : True,

        # Energy injection rate (m^2/s^3). amplitude is computed internally.
        # Do NOT set 'amplitude' when using epsilon_0.
        'epsilon_0'          : epsilon_0,

        # Forcing wavenumber band: forces L in (N_FORCING_MIN, N_FORCING_MAX)
        'n_total_forcing_min': N_FORCING_MIN,
        'n_total_forcing_max': N_FORCING_MAX,

        # Isotropic: no lat/lon mask, all m >= 1 forced.
        'do_localize'        : False,
        'zonal_forcing_min'  : 0,
    },

})

# -------------------------------------------------------------------------
# Run
# -------------------------------------------------------------------------
if __name__ == '__main__':

    cb.compile()

    # Spin-up: run for a long time to reach statistical equilibrium.
    # For barotropic turbulence, ~1000-2000 days is typically needed
    # before reliable statistics can be collected. Adjust as needed.
    exp.run(1, use_restart=False, num_cores=NCORES, overwrite_data=False)
    for i in range(2, 21):
       exp.run(i, num_cores=NCORES, overwrite_data=False)