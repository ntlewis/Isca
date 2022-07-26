import os

import numpy as np

from isca import SocratesCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE
from isca.util import exp_progress

NCORES = 12 
base_dir = os.path.dirname(os.path.realpath(__file__))
# a CodeBase can be a directory on the computer,
# useful for iterative development
cb = SocratesCodeBase.from_directory(GFDL_BASE)

# or it can point to a specific git repo and commit id.
# This method should ensure future, independent, reproducibility of results.
# cb = DryCodeBase.from_repo(repo='https://github.com/isca/isca', commit='isca1.1')

# compilation depends on computer specific settings.  The $GFDL_ENV
# environment variable is used to determine which `$GFDL_BASE/src/extra/env` file
# is used to load the correct compilers.  The env file is always loaded from
# $GFDL_BASE and not the checked out git repo.

# create an Experiment object to handle the configuration of model parameters
# and output diagnostics

exp = Experiment('soc_fixsst_annmean', codebase=cb)
exp.clear_rundir()

exp.inputfiles = [os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc'),
                  os.path.join(GFDL_BASE,'exp/test_cases/realistic_continents/input/era-spectral11_T63_96x192.out.nc'),
                  os.path.join(GFDL_BASE,'exp/test_cases/socrates_test/input/sst_annmean_amip.nc'), 
                  os.path.join(GFDL_BASE,'exp/test_cases/realistic_continents/input/siconc_clim_amip.nc')]

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')

#Write out diagnostics need for vertical interpolation post-processing
diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('dynamics', 'zsurf')

#Tell model which diagnostics to write
#diag.add_field('atmosphere', 'precipitation', time_avg=True)
#diag.add_field('atmosphere', 'rh', time_avg=True)
#diag.add_field('mixed_layer', 't_surf', time_avg=True)
#diag.add_field('mixed_layer', 'flux_t', time_avg=True) #SH
#diag.add_field('mixed_layer', 'flux_lhe', time_avg=True) #LH
#diag.add_field('dynamics', 'sphum', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
#diag.add_field('dynamics', 'vcomp', time_avg=True)
#diag.add_field('dynamics', 'omega', time_avg=True)
diag.add_field('dynamics', 'temp', time_avg=True)
#diag.add_field('dynamics', 'height', time_avg=True)
#diag.add_field('hs_forcing', 'local_heating', time_avg=True)


# diag.add_file('atmos_daily', 1, 'days', time_units='days')

# #Write out diagnostics need for vertical interpolation post-processing
# diag.add_field('dynamics', 'ps', time_avg=False)
# diag.add_field('dynamics', 'bk')
# diag.add_field('dynamics', 'pk')
# diag.add_field('dynamics', 'zsurf')

# #Tell model which diagnostics to write
# diag.add_field('atmosphere', 'precipitation', time_avg=False)
# #diag.add_field('atmosphere', 'rh', time_avg=True)
# diag.add_field('mixed_layer', 't_surf', time_avg=False)
# #diag.add_field('mixed_layer', 'flux_t', time_avg=True) #SH
# #diag.add_field('mixed_layer', 'flux_lhe', time_avg=True) #LH
# diag.add_field('dynamics', 'sphum', time_avg=False)
# diag.add_field('dynamics', 'ucomp', time_avg=False)
# diag.add_field('dynamics', 'vcomp', time_avg=False)
# diag.add_field('dynamics', 'omega', time_avg=False)
# diag.add_field('dynamics', 'temp', time_avg=False)
# diag.add_field('dynamics', 'height', time_avg=False)
# diag.add_field('hs_forcing', 'local_heating', time_avg=False)

#temperature tendency - units are K/s
#diag.add_field('socrates', 'soc_tdt_lw', time_avg=True) # net flux lw 3d (up - down)
#diag.add_field('socrates', 'soc_tdt_sw', time_avg=True)
#diag.add_field('socrates', 'soc_tdt_rad', time_avg=True) #sum of the sw and lw heating rates

#net (up) and down surface fluxes
#diag.add_field('socrates', 'soc_surf_flux_lw', time_avg=True)
#diag.add_field('socrates', 'soc_surf_flux_sw', time_avg=True)
#diag.add_field('socrates', 'soc_surf_flux_lw_down', time_avg=True)
#diag.add_field('socrates', 'soc_surf_flux_sw_down', time_avg=True)
#net (up) TOA and downard fluxes
# diag.add_field('socrates', 'soc_olr', time_avg=True)
# diag.add_field('socrates', 'soc_toa_sw', time_avg=True) 
# diag.add_field('socrates', 'soc_toa_sw_down', time_avg=True)

exp.diag_table = diag

#Define values for the 'core' namelist
exp.namelist = namelist = Namelist({
    'main_nml':{
     'days'   : 360,
     'hours'  : 0,
     'minutes': 0,
     'seconds': 0,
     'dt_atmos':600,
     'current_date' : [1,1,1,0,0,0],
     'calendar' : 'thirty_day'
    },
    'socrates_rad_nml': {
        'stellar_constant':1365.,
        'lw_spectral_filename':os.path.join(GFDL_BASE,'src/atmos_param/socrates/src/trunk/data/spectra/ga7/sp_lw_ga7'),
        'sw_spectral_filename':os.path.join(GFDL_BASE,'src/atmos_param/socrates/src/trunk/data/spectra/ga7/sp_sw_ga7'),
        'do_read_ozone': True,
        'ozone_file_name':'ozone_1990',
        'ozone_field_name':'ozone_1990',
        'dt_rad':3600,
        'store_intermediate_rad':True,
        'chunk_size': 16,
        'use_pressure_interp_for_half_levels':False,
        'tidally_locked':False,
        'solday':90
    }, 
    'idealized_moist_phys_nml': {
        'do_damping': True,
        'turb':True,
        'mixed_layer_bc':True,
        'do_virtual' :True,
        'do_simple': False,
        'roughness_mom':5.0e-03,
        'roughness_heat':1.0e-05,
        'roughness_moist':1.0e-05,            
        'two_stream_gray': False,     #Use the grey radiation scheme
        'do_socrates_radiation': True,
        'convection_scheme': 'SIMPLE_BETTS_MILLER', #Use simple Betts miller convection            
        'do_cloud_simple': False, 
        'land_option' : 'input',
        'land_file_name' : 'INPUT/era-spectral11_T63_96x192.out.nc',
        'land_roughness_prefactor' :10.0, 
        'do_local_heating':False, 
    },


    'vert_turb_driver_nml': {
        'do_mellor_yamada': False,     # default: True
        'do_diffusivity': True,        # default: False
        'do_simple': False,            # default: False
        'constant_gust': 1.0,          # default: 1.0
        'use_tau': False
    },
    
    'hs_forcing_nml': {
	'local_heating_option':'Polar', 
        'local_heating_srfamp':2.0, 
        'local_heating_ywidth':15.0, 
	'local_heating_ycenter':90.0, 
	'local_heating_vert_decay':5.0e3, 
	'local_heating_delta_p0':8.0e4, 
    }, 

    'diffusivity_nml': {
        'do_entrain':False,
        'do_simple': False,
    },

    'surface_flux_nml': {
        'use_virtual_temp': True,
        'do_simple': False,
        'old_dtaudv': False,    
	    'land_humidity_prefactor': 0.7,
    },

    'atmosphere_nml': {
        'idealized_moist_model': True
    },

    #Use a large mixed-layer depth, and the Albedo of the CTRL case in Jucker & Gerber, 2017
    'mixed_layer_nml': {
        'tconst' : 285.,
        'prescribe_initial_dist':True,
        'evaporation':True,  
        'land_option': 'input',              #Tell mixed layer to get land mask from input file
        'land_h_capacity_prefactor': 0.1,    #What factor to multiply mixed-layer depth by over land. 
        'albedo_value': 0.25,                #albedo value
        'land_albedo_prefactor': 1.3,        #What factor to multiply ocean albedo by over land     
        'do_qflux' : False, #Don't use the prescribed analytical formula for q-fluxes
        'do_read_sst' : True, #Read in sst values from input file
        'do_sc_sst' : True, #Do specified ssts (need both to be true)
        'sst_file' : 'sst_annmean_amip', #Set name of sst input file
        'specify_sst_over_ocean_only' : True, #Make sure sst only specified in regions of ocean.              
    },

    'qe_moist_convection_nml': {
        'rhbm':0.7,
        'tau_bm':7200.,
        'Tmin':120.,
        'Tmax':360.   
    },
    
    'lscale_cond_nml': {
        #'do_simple':True,
        'do_evap':True
    },
    
    'sat_vapor_pres_nml': {
           'do_simple':True,
           'construct_table_wrt_liq_and_ice':True
       },
    
    'damping_driver_nml': {
        'do_rayleigh': True,
        'trayfric': -0.5,              # neg. value: time in *days*
        'sponge_pbottom':  500., #Setting the lower pressure boundary for the model sponge layer in Pa.
        'do_conserve_energy': True,      
    },

    # FMS Framework configuration
    'diag_manager_nml': {
        'mix_snapshot_average_fields': False  # time avg fields are labelled with time in middle of window
    },

    'fms_nml': {
        'domains_stack_size': 600000                        # default: 0
    },

    'fms_io_nml': {
        'threading_write': 'single',                         # default: multi
        'fileset_write': 'single',                           # default: multi
    },

    'spectral_dynamics_nml': {
        'damping_order': 4,             
        'water_correction_limit': 200.e2,
        'reference_sea_level_press':1.0e5,
        'num_levels':48,      #How many model pressure levels to use
        'valid_range_t':[100.,800.],
        'initial_sphum':[2.e-6],
        'vert_coord_option':'input',
        'robert_coeff':0.03,
        'ocean_topog_smoothing': 0.0, 
        'num_fourier':63,
        'num_spherical':64,
        'lon_max':192, 
        'lat_max':96, 
    },

    'spectral_init_cond_nml':{
        'topog_file_name': 'era-spectral11_T63_96x192.out.nc',
        'topography_option': 'input'
    },
    
    'vert_coordinate_nml':{
        'bk': [0.00000,       0.00000,       0.00000,
               0.00000,       0.00000,       0.00000,
               0.00000,       0.00000,       0.00000,
               0.00000,       0.00000,       0.00000,
               0.00000,       0.00000,       0.00000,
               0.00000,       0.00000,       0.00000,
               0.00000,       0.00000,       0.00000,
               0.00000,       0.00000,       0.00000,
               0.00000,       0.00000,       0.00000,
               0.00000,       0.00000,       0.00000,
               0.00000,       0.00813,       0.03224,
               0.07128,       0.12445,       0.19063,
               0.26929,       0.35799,       0.45438,
               0.55263,       0.64304,       0.71703,
               0.77754,       0.82827,       0.87352,
               0.91502,       0.95235,       0.98511,
               1.00000  ],
        'pk': [40.00000,     100.00000,     200.00000,
              350.00000,     550.00000,     800.00000,
             1085.00000,    1390.00000,    1720.00000,
             2080.00000,    2470.00000,    2895.00000,
             3365.00000,    3890.00000,    4475.00000,
             5120.00000,    5830.00000,    6608.00000,
             7461.00000,    8395.00000,    9424.46289,
            10574.46900,   11864.80330,   13312.58850,
            14937.03770,   16759.70760,   18804.78670,
            21099.41250,   23674.03720,   26562.82650,
            29804.11680,   32627.31601,   34245.89759,
            34722.29104,   34155.20062,   32636.50533,
            30241.08406,   27101.45052,   23362.20912,
            19317.04955,   15446.17194,   12197.45091,
             9496.39912,    7205.66920,    5144.64339,
             3240.79521,    1518.62245,       0.00000,
                0.00000]
    },
    

})

#Lets do a run!
if __name__=="__main__":

        #cb.compile(debug=False)
        #Set up the experiment object, with the first argument being the experiment name.
        #This will be the name of the folder that the data will appear in.

        overwrite=False

        # spin up
        run_exp = exp.derive('soc_fixsst_annmean_T63')
        run_exp.run(1, use_restart=False, num_cores=NCORES, overwrite_data=overwrite)#, run_idb=True)
        for i in range(2,101):
#             if i == 41:
#                 diag = DiagTable()
#                 diag.add_file('atmos_6hrly', 6, 'hours', time_units='days')
                
#                 #Write out diagnostics need for vertical interpolation post-processing
#                 diag.add_field('dynamics', 'ps', time_avg=True)
#                 diag.add_field('dynamics', 'bk')
#                 diag.add_field('dynamics', 'pk')
#                 diag.add_field('dynamics', 'zsurf')

#                 #Tell model which diagnostics to write
#                 diag.add_field('mixed_layer', 't_surf', time_avg=True)
#                 diag.add_field('dynamics', 'sphum', time_avg=True)
#                 diag.add_field('dynamics', 'ucomp', time_avg=True)
#                 diag.add_field('dynamics', 'vcomp', time_avg=True)
#                 diag.add_field('dynamics', 'omega', time_avg=True)
#                 diag.add_field('dynamics', 'temp', time_avg=True)
#                 diag.add_field('dynamics', 'height', time_avg=True)
                
#                 run_exp.diag_table = diag 
                
            run_exp.run(i, num_cores=NCORES, overwrite_data=overwrite)
            
        
        
        
        #run_exp2 = exp.derive('soc_fixsst_annmean_noheat_daily50')
        #run_exp2.run(1, use_restart=True, num_cores=NCORES, overwrite_data=overwrite, 
        #             restart_file='/scratch/ntl203/isca_data/soc_fixsst_annmean_noheat/restarts/res0050.tar.gz')
        #for i in range(2,101):
#             if i == 51:
#                 diag = DiagTable()
#                 diag.add_file('atmos_6hrly', 6, 'hours', time_units='days')
                
#                 #Write out diagnostics need for vertical interpolation post-processing
#                 diag.add_field('dynamics', 'ps', time_avg=True)
#                 diag.add_field('dynamics', 'bk')
#                 diag.add_field('dynamics', 'pk')
#                 diag.add_field('dynamics', 'zsurf')

#                 #Tell model which diagnostics to write
#                 diag.add_field('mixed_layer', 't_surf', time_avg=True)
#                 diag.add_field('dynamics', 'sphum', time_avg=True)
#                 diag.add_field('dynamics', 'ucomp', time_avg=True)
#                 diag.add_field('dynamics', 'vcomp', time_avg=True)
#                 diag.add_field('dynamics', 'omega', time_avg=True)
#                 diag.add_field('dynamics', 'temp', time_avg=True)
#                 diag.add_field('dynamics', 'height', time_avg=True)
                
#                 run_exp2.diag_table = diag 
                
        #    run_exp2.run(i, num_cores=NCORES, overwrite_data=overwrite)
            
            
#        run_exp3 = exp.derive('soc_fixsst_annmean_localheat_daily50')
#        run_exp3.namelist['idealized_moist_phys_nml']['do_local_heating'] = True
#        
#        run_exp3.run(1, use_restart=True, num_cores=NCORES, overwrite_data=overwrite, 
#                     restart_file='/scratch/ntl203/isca_data/soc_fixsst_annmean_localheat/restarts/res0050.tar.gz')
#        for i in range(2,101):
#             if i == 51:
#                 diag = DiagTable()
#                 diag.add_file('atmos_6hrly', 6, 'hours', time_units='days')
                
#                 #Write out diagnostics need for vertical interpolation post-processing
#                 diag.add_field('dynamics', 'ps', time_avg=True)
#                 diag.add_field('dynamics', 'bk')
#                 diag.add_field('dynamics', 'pk')
#                 diag.add_field('dynamics', 'zsurf')

#                 #Tell model which diagnostics to write
#                 diag.add_field('mixed_layer', 't_surf', time_avg=True)
#                 diag.add_field('dynamics', 'sphum', time_avg=True)
#                 diag.add_field('dynamics', 'ucomp', time_avg=True)
#                 diag.add_field('dynamics', 'vcomp', time_avg=True)
#                 diag.add_field('dynamics', 'omega', time_avg=True)
#                 diag.add_field('dynamics', 'temp', time_avg=True)
#                 diag.add_field('dynamics', 'height', time_avg=True)
                
#                 run_exp3.diag_table = diag 
                
#            run_exp3.run(i, num_cores=NCORES, overwrite_data=overwrite)
