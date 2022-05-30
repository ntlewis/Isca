import numpy as np
from isca import IscaCodeBase, Experiment, DiagTable, Namelist, FailedRunError, log, GFDL_BASE
from math import pi




NDAYS = 20
NRUNS = 4  # total sim time = NDAYS*NRUNS
RESOLUTION = 'T42', 25
NCORES = 8
PERIOD = 10 * 86400 
OMEGA = 2 * 3.1415926 / PERIOD 
print(2 * pi / OMEGA, 86400.)


cb = IscaCodeBase.from_directory(GFDL_BASE)
cb.compile()
exp = Experiment('exp301_tl_frierson', codebase=cb)

namelist = Namelist({
    # model settings
    'main_nml': {
        'dt_atmos': 600,
        'days': NDAYS,
        'calendar': 'no_calendar',
    },
    'spectral_dynamics_nml': {
        # set level 7 at ~ 100hPa
        'damping_order': 2,
        'reference_sea_level_press': 1.0e5,  # default: 101325
        'valid_range_t': [100., 800.],  # default: (100, 500)
        'vert_coord_option': 'uneven_sigma',  # default: 'even_sigma'
        'scale_heights': 6.0,
        'exponent': 7.5,
        'surf_res': 0.5,
        'initial_sphum': [2.e-6],
        #        'p_sigma': .69,
        #        'p_press': .4
    },
    'astronomy_nml': {
        'ecc': 0.0,
        'obliq': 0.0
    },
    'idealized_moist_phys_nml': {
        'two_stream_gray': True,
        'do_rrtm_radiation': False,
        'convection_scheme': 'simple_betts_miller',
        'do_damping': True,
        'turb': True,
        #'lwet_convection': False,
        #'do_bm': True,
        'mixed_layer_bc': True,
        'do_virtual': False,
        'do_simple': True,
        # Roughness Lengths for Monin-Obukhov theory:
        # Baseline 3.21e-5 m
        # Ref:  Heng et al: Mon. Not. R. Astron. Soc [418] (2011)
        #       Frierson et al: J Atmos. Sci [63] (2006)
        # roughness_mom:
        #   Open water: 1e-4m
        #   Urban terrain: 1m
        'roughness_mom': 3.21e-05,  # default: 0.05
        'roughness_heat': 3.21e-05,  # default: 0.05
        'roughness_moist': 3.21e-05  # default: 0.05
    },
    'vert_turb_driver_nml': {
        'do_mellor_yamada': False,  # default: True
        'do_diffusivity': True,  # default: False
        'do_simple': True,  # default: False
        'constant_gust': 0.0,  # default: 1.0
        'use_tau': False
    },
    'damping_driver_nml': {
        'do_rayleigh': True,
        'do_conserve_energy': False,
        'trayfric': -3,  # neg. number = days
    },
    'diffusivity_nml': {
        'do_entrain': False,
        'do_simple': True,
    },
    'surface_flux_nml': {
        'use_virtual_temp': False,
        'do_simple': True,
        'old_dtaudv': True
    },
    'mixed_layer_nml': {
        'albedo_value': 0.31,  # for BOG, too much forcing so turn down temp
        'depth': 10.0,
        'tconst': 300.,  # default: 305.0
        'prescribe_initial_dist': True,
        'do_qflux': False,
    },
    'atmosphere_nml': {
        'idealized_moist_model': True
    },
    'spectral_init_cond_nml': {
        'initial_temperature': 300.0
    },
    'dry_convection_nml': {
        'gamma': 0.7,
        'tau': 7200.0 * 2,
    },
    'constants_nml': {
        'omega': OMEGA,
        'orbital_period': PERIOD,  # tidally locked
    },
    'astronomy_nml': {
        'ecc': 0.0,
        'obliq': .0,  # planetary obliquity
        'per': 0.0  # begin at autumn equinox
    },
    'two_stream_gray_rad_nml': {
        'rad_scheme': 'frierson',
        'do_seasonal': True,
        'atm_abs': 0.0,  # no SW absorption in the atmosphere
        'linear_tau': 0.0,  # no LW absorption in the stratosphere
        'ir_tau_pole': 6.0,
        'ir_tau_eq': 6.0,
        'wv_exponent': 1.0,
    },
    'qe_moist_convection_nml': {
        'rhbm': 0.7,
        'Tmin': 130.,
        'Tmax': 350.
    },
    # framework and IO config
    'diag_manager_nml': {
        # diagmanager gives a warning if you don't set this
        'mix_snapshot_average_fields': False
    },
    'fms_nml': {
        'domains_stack_size': 600000  # default: 0
    },
    'fms_io_nml': {
        'threading_write': 'single',  # default: multi
        'fileset_write': 'single',  # default: multi
    },
})


diag = DiagTable()
diag.add_file('daily', 1, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('atmosphere', 'precipitation', time_avg=True)
diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('dynamics', 'sphum', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'vor', time_avg=True)
diag.add_field('dynamics', 'div', time_avg=True)

diag.add_field('two_stream', 'coszen', time_avg=True)
diag.add_field('two_stream', 'olr', time_avg=True)

exp.namelist = namelist
exp.diag_table = diag
exp.set_resolution(*RESOLUTION)

exp.run(1, use_restart=False, num_cores=NCORES)
