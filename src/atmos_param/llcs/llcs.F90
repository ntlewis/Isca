! The Lambert-Lewis Convection Scheme (LLCS)
!
! Moist convective adjustment scheme with a dry limit. Convection can begin at any
! model levels, and multiple isolated convective 'events' can occur in a single
! column. 
!
! The scheme is documented in Appendix A of Lambert et al., 2020. Continuous
! Structural Parametrization: A method for representing different model
! parameterizations within one structure demonstrated for atmospheric convection,
! J. Adv. Model Earth Sys., in review.
!
! Code author: Neil Lewis (Oxford), with contributions from: Ian Boutle (UK Met Office)
! and Denis Sergeev (Exeter). 

MODULE llcs_mod

! Description:
!  Module containing all subroutines for the Lambert-Lewis
!  Convection Scheme (LLCS)
!  Further description of the scheme can be found below



USE              fms_mod, ONLY: file_exist, error_mesg, open_file, &
                              check_nml_error, mpp_pe, FATAL,    & 
                              NOTE, close_file, uppercase 
USE   sat_vapor_pres_mod, ONLY: escomp              
USE        constants_mod, ONLY: cp_air, grav, rdgas, rvgas, kappa, hlv, cp_vapor


IMPLICIT NONE
PRIVATE


PUBLIC llcs_init, llcs

! version number 
character(len=128) :: version = '$Id: llcs.F90,v 1.0 2020/06/12'
character(len=128) :: tag = '$Name:  $'
character(len=4) :: mod_name = 'llcs'

! local initialisation flag
LOGICAL :: do_init = .TRUE. 

! -- OPTIONS -- !
INTEGER, PARAMETER :: llcs_opt_all_rain = 0 !NTL: different cloud options for LLCS
                                            ! default ALL_RAIN, i.e. there are no clouds and all 
                                            ! condensate removed as precipitation 
INTEGER, PARAMETER :: llcs_opt_all_cloud = 1   !NTL: different cloud options for LLCS                                     
INTEGER, PARAMETER :: llcs_opt_crit_condens = 2   !NTL: different cloud options for LLCS 
INTEGER, PARAMETER :: llcs_opt_const_frac = 3   !NTL: different cloud options for LLCS 
character(len=256) :: cloud_option = 'ALL_RAIN' !overwritten in namelist 
INTEGER :: llcs_cloud_precip = 0 ! overwritten in llcs_init 


! -- PARAMETERS -- !
! Parameters (of the scheme / model) overwritten in namelist 
REAL :: llcs_detrain_coef = 0.6 !NTL: what is this? something to do with clouds 
REAL :: llcs_rhcrit = 0.8 ! critical relative humidity for triggering moist convection 
REAL :: llcs_timescale = 3600.0 ! timescale for moist convective adjustment 
REAL :: llcs_rain_frac = 0.5 !NTL: what is this? something to do with clouds 
REAL :: fac_qsat = 0.500 !NTL: what is this? something to do with clouds 
REAL :: qlmin = 3.0000e-4 !NTL: what is this? something to do with clouds 
REAL :: mparwtr = 1.5000e-3    !NTL: what is this? something to do with clouds 

! option for qsat_calc and rsat_calc as for rh_calc in idealized moist phys 
LOGICAL :: do_simple = .FALSE. !                             

! other parameters (not overwritten in namelist )
REAL, PARAMETER :: p_zero = 1.e5 ! reference pressure for conversion between theta and temperature 
REAl :: repsilon = 0.0 ! rdgas/rvgas, initialised in llcs_init
REAL :: one_minus_epsilon = 0.0 ! 1 - repsilon, initialised in llcs_init
REAL :: kappa_use = 0.0 ! modified kappa (kappa_use = kappa * alpha) 
                        ! initialised in llcs_init
REAL :: alpha = 1.0 ! DONT CHANGE FROM 1.0 unless the atmosphere is dry 


namelist /llcs_nml/ llcs_detrain_coef, llcs_rhcrit, llcs_timescale, llcs_rain_frac, & 
                    fac_qsat, qlmin, mparwtr, cloud_option, do_simple, alpha 

CONTAINS



SUBROUTINE llcs_init ()

    !-----------------------------------------------------------------------
    !
    !        initialization of LLCS convection scheme 
    !
    !-----------------------------------------------------------------------
  
    INTEGER ::  unit, io, ierr

    !----------- read namelist ---------------------------------------------

    if (file_exist('input.nml')) then
       unit = open_file (file='input.nml', action='read')
       ierr = 1
       do while (ierr /= 0)
          read  (unit, nml=llcs_nml, iostat=io, end=10)
          ierr = check_nml_error(io, 'llcs_nml')
       end do
10     call close_file(unit)
    endif

    !---------- output namelist --------------------------------------------

    unit = open_file (file='logfile.out', action='append')
    if ( mpp_pe() == 0 ) then
       write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
       write (unit,nml=llcs_nml)
    endif
    call close_file(unit)

    do_init = .FALSE.
    
    ! initialise repsilon and one_minus_epsilon
    repsilon = rdgas / rvgas 
    one_minus_epsilon = 1 - repsilon

    ! initialise kappa_use 
    kappa_use = kappa * alpha 

    ! select llcs cloud option 
    if(uppercase(trim(cloud_option)) == 'ALL_RAIN') then
      llcs_cloud_precip = llcs_opt_all_rain
      call error_mesg('llcs','all_rain mode selected for LLCS', NOTE)
    elseif(uppercase(trim(cloud_option)) == 'ALL_CLOUD') then
      llcs_cloud_precip = llcs_opt_all_cloud
      call error_mesg('llcs','all_cloud mode selected for LLCS.'//                 & 
                             'This is NOT a valid option, as this means LLCS'//    &
                             'will pass some condensate to cloud liquid water, '// &
                             'which is currently not implemented in Isca. FATAL.', &
                             FATAL)
    elseif(uppercase(trim(cloud_option)) == 'CRIT_CONDENS') then
      llcs_cloud_precip = llcs_opt_crit_condens
      call error_mesg('llcs','crit_condens mode selected for LLCS.'//              & 
                             'This is NOT a valid option, as this means LLCS'//    &
                             'will pass some condensate to cloud liquid water, '// &
                             'which is currently not implemented in Isca. FATAL.', &
                             FATAL)
    elseif(uppercase(trim(cloud_option)) == 'CONST_FRAC') then
      llcs_cloud_precip = llcs_opt_const_frac
      call error_mesg('llcs','const_frac mode selected for LLCS.'//                & 
                             'This is NOT a valid option, as this means LLCS'//    &
                             'will pass some condensate to cloud liquid water, '// &
                             'which is currently not implemented in Isca. FATAL.', &
                             FATAL)
    else 
      call error_mesg('llcs','No cloud option selected for LLCS. Defaulting to ALL_RAIN.', NOTE)
    end if 



       
END SUBROUTINE llcs_init




! ----------------------------------------------------------------------
! -- ISCA WRAPPER -- !
! -- Call to convection scheme is CALL convection_scheme(....)
! ----------------------------------------------------------------------

SUBROUTINE llcs(                                                       &
    ! Inputs
    temp_s_arr, q_s_arr, qcl_s_arr, p_half_arr, p_full_arr, timestep,           &
    ! Outputs
    temp_f_arr, q_f_arr, qcl_inc_arr, cfl_arr, rain_arr)



IMPLICIT NONE

! Inputs
REAL, INTENT(IN) :: temp_s_arr   (:,:,:)
REAL, INTENT(IN) :: q_s_arr   (:,:,:)
REAL, INTENT(IN) :: qcl_s_arr (:,:,:)
REAL, INTENT(IN) :: p_full_arr(:,:,:)
REAL, INTENT(IN) :: p_half_arr(:,:,:)
REAL, INTENT(IN) :: timestep

! Outputs
REAL, INTENT(OUT) :: temp_f_arr    (:,:,:)
REAL, INTENT(OUT) :: q_f_arr    (:,:,:)
REAL, INTENT(OUT) :: qcl_inc_arr(:,:,:)
REAL, INTENT(OUT) :: cfl_arr    (:,:,:)
REAL, INTENT(OUT) :: rain_arr   (:,:)

! Locals
INTEGER :: i
INTEGER :: j
INTEGER :: nlevels
integer :: rows 
integer :: row_length 



rows = SIZE(temp_s_arr, 2)
row_length = SIZE(temp_s_arr, 1)
nlevels = SIZE(temp_s_arr, 3)


! Check whether initialization has been completed
    if (do_init) call error_mesg ('llcs',  &
         'llcs_init has not been called.', FATAL)

DO j = 1, rows
  DO i = 1, row_length
    CALL convection_scheme(                                                       &
        ! Inputs
        nlevels, temp_s_arr(i, j, :), q_s_arr(i, j, :), qcl_s_arr(i, j, :),       &
        p_half_arr(i, j, :), p_full_arr(i, j, :),  timestep,                      &
        ! Outputs
        temp_f_arr(i, j, :), q_f_arr(i, j, :), qcl_inc_arr(i, j, :),              &
        cfl_arr(i, j, :), rain_arr(i, j))
  END DO
END DO


END SUBROUTINE llcs

! ----------------------------------------------------------------------
! -- THE MAIN LLCS SUBROUTINE -- !
! ----------------------------------------------------------------------
SUBROUTINE convection_scheme(                                                  &
    ! Inputs
    nlevels, temp_start, q_start, qcl_start, p_half, p_full, timestep,         &
    ! Outputs
    temp_final, q_final, qcl_inc, cloud_frac, rain_sec)
!------------------------------------------------------------------------------
! Run Lambert-Lewis convection scheme.
!
! Required inputs are:
! potential temperature, specific humidity, pressure
! on full model levels, pressure on half model levels
! Potential temperature and specific humidity have (nlevels) levels.
! The arrays are set up such that, for example, temp_start(1) is
! temperature at the lowest model level and temp_start(nlevels) is
! temperature at the top of atmosphere.
!
! The convection scheme will output final temperature and
! final specific humidity, along with rainfall / second (kgm^-2s^-1)
! within the column.
!------------------------------------------------------------------------------
IMPLICIT NONE

! Inputs
INTEGER, INTENT(IN) :: nlevels ! number of full model levels
REAL, INTENT(IN) :: temp_start(:) ! inital temp. at full model levels.
REAL, INTENT(IN) :: q_start(:) ! initial specific humidity
REAL, INTENT(IN) :: qcl_start(:) ! initial cloud liquid content
REAL, INTENT(IN) :: p_full(:) ! pressure at full model levels
REAL, INTENT(IN) :: p_half(:) ! pressure at half model levels
REAL, INTENT(IN) :: timestep 
! Outputs
REAL, INTENT(OUT) :: temp_final(:) ! final temp. [K]
REAL, INTENT(OUT) :: q_final(:) ! final specific humidity [kg kg-1]
REAL, INTENT(OUT) :: qcl_inc(:) ! increment of cloud liquid water [kg kg-1]
REAL, INTENT(OUT) :: cloud_frac(:) ! cloud fraction [1]
REAL, INTENT(OUT) :: rain_sec ! convective precipitation [kg m-2 s-1]
! Local variables
INTEGER :: j_lev ! level loop index
INTEGER :: convbase ! the lowest level from which any form of mixing occurs.
INTEGER :: cloudbase ! level from which adjustment to a moist adiabat can occur
INTEGER :: final_level ! last level of a convective event.
                       ! Set to -1 whilst the event is in progress
REAL :: theta_start      (1:nlevels) ! theta calculated from initial temperature 
REAL :: theta_adjust     (1:nlevels) ! theta after correction to reference
                                     ! profile
REAL :: theta_noq_adjust (1:nlevels) ! dummy theta profile where only
                                     ! dry corr. happens, used for enthalpy
                                     ! conservation.
REAL :: theta_adj_rlx    (1:nlevels) ! theta_adjust corrected to account for
                                     ! timestep over timescale
REAL :: theta_noq_adj_rlx(1:nlevels) ! dry correction theta adjusted for
                                     ! timestep over timescale
REAL :: q_adjust         (1:nlevels) ! spec. hum. after correction to q_sat
REAL :: q_adj_rlx        (1:nlevels) ! q_adjust adjusted for timestep/timescale
REAL :: q_conserve       (1:nlevels) ! conserved sp. hum. in the column, other
                                     ! than moisture lost through rainfall
REAL :: q_dummy          (1:nlevels) ! dummy q profile to ensure moisture is
                                     ! not mixed twice
REAL :: delta_p          (1:nlevels) ! Change in pressure between the top and
                                     ! the bottom of a layer
REAL :: xmin             (1:nlevels) ! critical cloud condensate content
REAL :: rh_use ! RH adjustment limit, currently is constant and equal to 1.0

LOGICAL :: convection_flag ! indicator that convection takes place
LOGICAL :: moist_has_happened ! flag indicating whether moist convection occurs
LOGICAL :: moist_trigger ! flag for triggering moist convection

delta_p = p_half(1:nlevels) - p_half(2:nlevels+1)

! Initialise variables
CALL initialise(                                                               &
    ! Inputs
    nlevels, temp_start, q_start, p_full,                                      &
    ! Outputs
    final_level,                                                               &
    theta_start, theta_adjust, theta_noq_adjust,                               &
    theta_adj_rlx, theta_noq_adj_rlx, temp_final,                              &
    q_adjust, q_adj_rlx, q_conserve, q_dummy, q_final, qcl_inc,                &
    cloud_frac, xmin, rain_sec,                                                &
    moist_has_happened, convection_flag)

! Calculate critical cloud condensate
IF (llcs_cloud_precip == llcs_opt_crit_condens) THEN
  CALL critical_condensate(                                                    & 
      ! Inputs
      nlevels, theta_start, q_start, p_full,                                   &
      ! Outputs
      xmin)
END IF

! Diagnose and then perform convection.
DO j_lev = 1, nlevels-1
  IF (final_level == j_lev) THEN
    CALL diagnose(                                                             &
        ! Inputs
        nlevels, j_lev, theta_adjust, q_adjust, p_full,                        &
        ! In/out
        convection_flag,                                                       &
        ! Outputs
        convbase, cloudbase, final_level, moist_trigger, rh_use)
    IF (convbase /=  -1) THEN
      ! Convection: dry or moist
      IF (cloudbase == -1) THEN
        CALL dry_convection(                                                   &
            ! Inputs
            nlevels, convbase,                                                 &
            ! Inputs/outputs
            theta_adjust, theta_noq_adjust, q_adjust,                          &
            ! Outputs
            final_level)
      ELSE
        CALL moist_convection(                                                 &
            ! Inputs
            nlevels, convbase, cloudbase, rh_use, p_full,                      &
            ! Inputs/outputs
            theta_adjust, theta_noq_adjust, q_adjust, q_dummy,                 &
            moist_has_happened, moist_trigger,                                 &
            ! Outputs
            final_level)
      END IF
    END IF
  END IF
END DO

IF (convection_flag) THEN
  ! Convection has happened so we need to perform the following steps.

  ! Timescale relaxation
  ! In case of moist convection this step is performed twice:
  !   - for the actual theta profile,
  !   - for the 'dry convection only' profile used for enthalpy conservation
  CALL time_relax(                                                             &
      ! Inputs
      theta_start, theta_adjust, q_start, q_adjust, timestep,                  &
      ! Inputs/outputs
      theta_adj_rlx, q_adj_rlx)
  IF (moist_has_happened) THEN
    ! Second call for enthalpy conservation calculation
    CALL time_relax(                                                           &
        ! Inputs
        theta_start, theta_noq_adjust, q_start, q_adjust, timestep,            &
        ! Inputs/outputs
        theta_noq_adj_rlx, q_dummy)
  END IF

  ! Conserve Moisture
  CALL conserve_moisture(                                                      &
      ! Inputs
      nlevels, p_full, delta_p, q_start, q_adj_rlx,                            &
      ! Outputs
      q_conserve)

  IF (moist_has_happened) THEN
    ! Latent heating & equiv. amount of instantaneous rainfall
    CALL rain(                                                                 &
        ! Inputs
        nlevels, theta_start, theta_noq_adj_rlx, q_conserve, qcl_start, xmin,  &
        p_full, delta_p, timestep,                                             &
        ! Inputs/outputs
        theta_adj_rlx,                                                         &
        ! Outputs
        rain_sec, q_final, qcl_inc, cloud_frac)
  END IF

  ! Conserve enthalpy
  IF (moist_has_happened) THEN
    CALL conserve_enthalpy(                                                    &
        ! Inputs
        nlevels, theta_start, theta_adj_rlx, theta_noq_adj_rlx, p_full,        &
        delta_p,                                                               &
        ! Outputs
        temp_final)
  ELSE
    CALL conserve_enthalpy(                                                    &
        ! Inputs
        nlevels, theta_start, theta_adj_rlx, theta_adj_rlx, p_full, delta_p,   &
        ! Outputs
        temp_final)
  END IF

END IF

! Note that if no convection happens, then theta_final, q_final and
! rainsec are left as they are initialised in 'initialise', so in this
! case theta_final = theta_start, q_final = q_start and rainsec = 0.
END SUBROUTINE convection_scheme
! ----------------------------------------------------------------------
! -- END OF MAIN LLCS SUBROUTINE -- !
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! -- ROUTINES CALLED BY THE CONVECTION SCHEME -- !
! ----------------------------------------------------------------------
SUBROUTINE initialise(                                                         &
    ! Inputs
    nlevels, temp_start, q_start, p_full,                                      &
    ! Outputs
    final_level,                                                               &
    theta_start, theta_adjust, theta_noq_adjust,                               &
    theta_adj_rlx, theta_noq_adj_rlx, temp_final,                              &
    q_adjust, q_adj_rlx, q_conserve, q_dummy, q_final, qcl_inc,                &
    cloud_frac, xmin, rain_sec,                                                &
    moist_has_happened, convection_flag)
!------------------------------------------------------------------------------
! Initialise working arrays by setting them equal to input arrays.
!
! Further, the convbase and cloudbase flags are set to -1, as no
! convection has happened yet.
!
! Finally, rainsec is set to 0.
!------------------------------------------------------------------------------
IMPLICIT NONE

! Inputs
INTEGER, INTENT(IN) :: nlevels ! number of full model levels

REAL, INTENT(IN) :: temp_start(1:nlevels) ! inital temp. at full levels.
REAL, INTENT(IN) :: q_start(1:nlevels) ! initial specific humidity
REAL, INTENT(IN) :: p_full(1:nlevels) ! pressure on full model levels 

! Outputs
INTEGER, INTENT(OUT) :: final_level

REAL, INTENT(OUT) :: theta_start      (1:nlevels)
REAL, INTENT(OUT) :: theta_adjust     (1:nlevels)
REAL, INTENT(OUT) :: theta_noq_adjust (1:nlevels)
REAL, INTENT(OUT) :: theta_adj_rlx    (1:nlevels)
REAL, INTENT(OUT) :: theta_noq_adj_rlx(1:nlevels)
REAL, INTENT(OUT) :: temp_final      (1:nlevels)
REAL, INTENT(OUT) :: q_adjust         (1:nlevels)
REAL, INTENT(OUT) :: q_adj_rlx        (1:nlevels)
REAL, INTENT(OUT) :: q_conserve       (1:nlevels)
REAL, INTENT(OUT) :: q_dummy          (1:nlevels)
REAL, INTENT(OUT) :: q_final          (1:nlevels)
REAL, INTENT(OUT) :: qcl_inc          (1:nlevels)
REAL, INTENT(OUT) :: cloud_frac       (1:nlevels)
REAL, INTENT(OUT) :: xmin             (1:nlevels)
REAL, INTENT(OUT) :: rain_sec

LOGICAL, INTENT(OUT) :: convection_flag
LOGICAL, INTENT(OUT) :: moist_has_happened

call calc_theta(temp_start, p_full, theta_start) 
theta_adjust = theta_start
theta_noq_adjust = theta_start
theta_adj_rlx = theta_start
theta_noq_adj_rlx = theta_start
temp_final = temp_start

q_adjust = q_start
q_adj_rlx = q_start
q_conserve = q_start
q_final = q_start
q_dummy = q_start

qcl_inc = 0.0
cloud_frac = 0.0
xmin = 0.0
rain_sec = 0.0

final_level = 1
moist_has_happened = .FALSE.
convection_flag = .FALSE.

END SUBROUTINE initialise


SUBROUTINE diagnose(                                                           &
     ! Inputs
     nlevels, current_level, theta_start, q_start, p_full,                     &
     ! In/out
     convection_flag,                                                          &
     ! Outputs
     convbase, cloudbase, final_level, moist_trigger, rh_use)
!------------------------------------------------------------------------------
! Diagnose on which levels dry or moist convection occurs.
!
! This subroutine decides whether it is possible for convection to
! originate on the model level under consideration, and if so, whether
! at any point above this level moist convection is possible (i.e
! whether a parcel rising from this level would saturate).
!
! The criterion for whether mixing is possible is the dry instability:
! where d_theta / d_z < 0.
!
! Further the requirement for a rising parcel to saturate is
! q(convbase) > rh_crit * q_sat(level where it saturates).
!------------------------------------------------------------------------------

IMPLICIT NONE

! Inputs
INTEGER, INTENT(IN) :: nlevels
INTEGER, INTENT(IN) :: current_level

REAL, INTENT(IN) :: theta_start(1:nlevels)
REAL, INTENT(IN) :: q_start    (1:nlevels)
REAL, INTENT(IN) :: p_full     (1:nlevels)

! Input/Outputs
LOGICAL, INTENT(IN OUT) :: convection_flag

! Outputs
INTEGER, INTENT(OUT) :: convbase
INTEGER, INTENT(OUT) :: cloudbase
INTEGER, INTENT(OUT) :: final_level

REAL, INTENT(OUT) :: rh_use

LOGICAL, INTENT(OUT) :: moist_trigger

! Local variables
INTEGER :: j_lev
REAL :: temp
REAL :: q_sat
REAL :: q_crit
!REAL :: rh

! initialise cloudbase
cloudbase = -1
convbase = -1
final_level = -1
rh_use = 1.0
moist_trigger = .FALSE.

temp = theta_start(current_level) * (p_full(current_level)/p_zero) ** (kappa_use)

CALL qsat_calc(p_full(current_level), temp, q_sat)
q_crit = llcs_rhcrit * q_sat
moist_trigger = (q_start(current_level) >= q_crit)
IF ( (theta_start(current_level) >= theta_start(current_level+1)) .OR.         &
     moist_trigger ) THEN
  ! Thermal instability or moist instability

  IF (q_start(current_level) >= 0.0) THEN

    convbase = current_level

    ! Convection happens, change flag
    convection_flag = .TRUE. ! TODO: shouldn't it be outside this IF-statement?

    DO j_lev = current_level, nlevels - 1
      temp = theta_start(j_lev) * (p_full(j_lev)/p_zero) ** (kappa_use)

      CALL qsat_calc(p_full(current_level), temp, q_sat)

      q_crit = llcs_rhcrit * q_sat
      IF (q_start(j_lev) >= q_crit) THEN
        cloudbase = MAX(j_lev, convbase + 1) ! convbase cannot be cloudbase
        ! IF (cloudbase == convbase) THEN
        !   cloudbase = convbase + 1 ! convbase cannot be cloudbase
        ! END IF
        EXIT
      END IF
    END DO

    ! Allows reference profile RH to adopt RH of layer where moist
    ! convection begins.
    ! This is commented out for now pending further testing, but the code
    ! is left here for ease of future development
    ! IF (cloudbase /= -1) THEN
    !   temp = theta_start(cloudbase) * (p_full(cloudbase)/p_zero) ** (kappa_use)
    !   CALL qsat_calc(p_full(cloudbase), temp, q_sat)
    !   rh = q_start(cloudbase) / q_sat
    !   IF ((rh > rh_use) .AND. (rh <= 1.0)) THEN
    !     rh_use = rh
    !   ELSE IF (rh > 1.0) THEN
    !     rh_use = 1.0
    !   END IF
    ! END IF

  ELSE
    final_level = current_level + 1
  END IF
ELSE
  final_level = current_level + 1
END IF

END SUBROUTINE diagnose


SUBROUTINE dry_convection(                                                     &
    ! Inputs
    nlevels, convbase,                                                         &
    ! Inputs/outputs
    theta_adjust, theta_noq_adjust, q_adjust,                                  &
    ! Outputs
    final_level)
!------------------------------------------------------------------------------
! Perform dry convection.
!
! This routine increases the potential temperature on model levels to
! those given by the dry adiabat, where there is dry instability.
!
! We also perform the same adjustment for theta_noq_adjust, a dummy 'dry
! convection only' profile that we use for enthalpy conservation.
!
! Specific humidity is also mixed in a similar way, by setting
! q_adjust(level_under_consideration) to q_start(convbase).
!
! Both of these adjustments will then be reduced to account for
! timestep / timescale.
!
! Returns
! final_level: the level where the convective event ends.
!   It is from this level that the routine 'diagnose' will be re-called, in
!   an effort to find additional thermal instabilities that would permit
!   convection.
!------------------------------------------------------------------------------

IMPLICIT NONE

! Inputs
INTEGER, INTENT(IN) :: nlevels
INTEGER, INTENT(IN) :: convbase

! Inputs/outputs
REAL, INTENT(IN OUT) :: theta_adjust    (1:nlevels)
REAL, INTENT(IN OUT) :: theta_noq_adjust(1:nlevels)
REAL, INTENT(IN OUT) :: q_adjust        (1:nlevels)

! Outputs
INTEGER, INTENT(OUT) :: final_level

! Local variables
INTEGER :: j_lev

final_level = -1 ! if this remains unchanged then, in circumstances
                 ! where dry_convection has been called by moist_
                 ! convection, it signals to moist_convection that
                 ! the event is able to proceed. I.e the dry mixing
                 ! reached the level cloudbase.

DO j_lev = convbase + 1, nlevels
  IF ( (theta_adjust(j_lev-1) >= theta_adjust(j_lev)) .AND.                    &
       (q_adjust(j_lev-1) >=  0.0) ) THEN
    theta_adjust(j_lev) = theta_adjust(j_lev-1)
    theta_noq_adjust(j_lev) = theta_noq_adjust(j_lev-1)
    IF ( (q_adjust(j_lev-1) >= 0.0) .AND.                                      &
         (q_adjust(j_lev-1) > q_adjust(j_lev)) ) THEN
      q_adjust(j_lev) = q_adjust(j_lev-1)
    END IF
  ELSE
    final_level = j_lev
    EXIT
  END IF
END DO

END SUBROUTINE dry_convection

! ----------------------------------------------------------------------

SUBROUTINE moist_convection(                                                   &
    ! Inputs
    nlevels, convbase, cloudbase, rh_use, p_full,                              &
    ! Inputs/outputs
    theta_adjust, theta_noq_adjust, q_adjust, q_dummy,                         &
    moist_has_happened, moist_trigger,                                         &
    ! Outputs
    final_level)
!------------------------------------------------------------------------------
! Perform moist convection.
!
! This routine is called if we know that a level cloudbase exists
! such that a parcel of air rising from convbase would saturate.
!
! In this situation the theta profile is adjustted beyond that of a
! dry reference profile to the reference profile for the moist
! pseudo-adiabat. This is as a result of the latent heating provided
! by saturation.
!
! The specific humidity within the columns is mixed as before, unless
! this would lead to a level having specific humidity greater than
! saturation specific humidity, in which case we set specific humidity
! equal to saturation specific humidity.
!
! If cloudbase-convbase >= 2 then we perform dry convection until
! we reach the level cloudbase.
!
! If dry convection ceases to be possible before we reach the level
! cloudbase then we claim that a parcel rising from convbase does not
! reach a level at which it saturates, and so we end the convective event
! accordingly without actually performing moist convection. Moist
! covection may still be possible in such a scenario if, when we re-call
! diagnose from 'final_level', it finds another level convbase with
! accompanying level cloudbase.
!
! Finally it returns an integer, moist_has_happened, which signals whether
! moist convection has infact taken place. This is initalised in the
! routine 'initialise' to be moist_has_happened = -1. If at any point
! moist convection takes place this is changed to 1.
!------------------------------------------------------------------------------
IMPLICIT NONE

! Inputs
INTEGER, INTENT(IN) :: nlevels
INTEGER, INTENT(IN) :: convbase
INTEGER, INTENT(IN) :: cloudbase

REAL, INTENT(IN) :: rh_use
REAL, INTENT(IN) :: p_full(1:nlevels)

! Input/outputs
REAL, INTENT(IN OUT) :: theta_adjust    (1:nlevels)
REAL, INTENT(IN OUT) :: theta_noq_adjust(1:nlevels)
REAL, INTENT(IN OUT) :: q_adjust        (1:nlevels)
REAL, INTENT(IN OUT) :: q_dummy         (1:nlevels)

LOGICAL, INTENT(IN OUT) :: moist_has_happened
LOGICAL, INTENT(IN OUT) :: moist_trigger

! Outputs
INTEGER, INTENT(OUT) :: final_level

! Locals
INTEGER :: j_lev
REAL :: temp(1:nlevels)
REAL :: q_sat
REAL :: r_sat
REAL :: temp_v
REAL :: gamma_s
REAL :: gamma_p
REAL :: rho
REAL :: q_adjust_save

final_level = -1 ! final_level = -1 allows moist convection to happen.

CALL calc_temp(theta_adjust, p_full, temp)

IF ( (cloudbase-convbase) >= 2 ) THEN
  ! Note that in this situation we only want to perform dry convection
  ! up until the level below cloudbase, so we use cloudbase-1 as the
  ! 'nlevels' input.
  CALL dry_convection(                                                         &
      ! Inputs
      cloudbase-1, convbase,                                                   &
      ! Inputs/outputs
      theta_adjust, theta_noq_adjust, q_adjust,                                &
      ! Outputs
      final_level)
  CALL calc_temp(theta_adjust, p_full, temp)
END IF

IF (final_level == -1) THEN
  ! Dry convection reached level below cloudbase, so we do moist convection.
  DO j_lev = cloudbase, nlevels
    IF ( (theta_adjust(j_lev-1) >= theta_adjust(j_lev))                        &
        .OR. moist_trigger ) THEN

      moist_trigger = .FALSE. ! reset moist trigger

      IF (q_adjust(j_lev-1) >= 0.0) THEN
        ! Adjust specific humidity    CHANGES HERE REMOVE RH_CRITS
        CALL qsat_calc(p_full(j_lev), temp(j_lev), q_sat)
        q_sat = rh_use * q_sat

        ! Save pre-adjusted spec. hum. as local scalar
        q_adjust_save = q_adjust(j_lev)

        ! Adjust spec. hum. at the current level
        q_adjust(j_lev) = MAX(q_adjust(j_lev), q_sat)

        ! Next, adjust potential temperature

        CALL rsat_calc(p_full(j_lev), temp(j_lev), r_sat)
        r_sat = rh_use * r_sat

        ! Saturation dT / dz
        gamma_s = -1 * grav *                                                         &
            ( (1+r_sat) * (1+(hlv*r_sat) / (rdgas * temp(j_lev-1))) )               &
            / (cp_air + r_sat*cp_vapor +                                              &
            (repsilon + r_sat) * ((r_sat*hlv**2)/(rdgas*temp(j_lev-1)**2)))

        CALL rsat_calc(p_full(j_lev-1), temp(j_lev-1), r_sat)
        r_sat = rh_use * r_sat

        temp_v = temp(j_lev-1) * (1 + r_sat / repsilon) / (1 + r_sat)
        ! Virtual temperature temp_v used for density of moist air
        rho = p_full(j_lev-1) / (rdgas * temp_v) ! moist air density
        gamma_p = (1.0 / (rho * grav)) * gamma_s ! saturated dT/dp
        IF ( (temp(j_lev-1) +                                                  &
             (gamma_p * (p_full(j_lev-1)-p_full(j_lev)))) <= temp(j_lev) ) THEN
          final_level = j_lev
          q_adjust(j_lev) = q_adjust_save
          EXIT
        END IF

        ! Reference profile temperature
        temp(j_lev) = temp(j_lev-1) +                                          &
            gamma_p * (p_full(j_lev-1) - p_full(j_lev))
        ! adjust potential temperature accordingly
        theta_adjust(j_lev) = temp(j_lev) * (p_zero / p_full(j_lev))**kappa_use

        moist_has_happened = .TRUE.

        ! Adjust dummy dry potential temperature profile
        theta_noq_adjust(j_lev) = MAX(theta_noq_adjust(j_lev-1),               &
                                      theta_noq_adjust(j_lev))
      ELSE
        final_level = j_lev
        EXIT
      END IF

    ELSE ! record final level
      final_level = j_lev
      EXIT
    END IF
  END DO
END IF

END SUBROUTINE moist_convection


SUBROUTINE time_relax(                                                         &
    ! Inputs
    theta_start, theta_adjust, q_start, q_adjust, timestep,                     &
    ! Inputs/outputs
    theta_adj_rlx, q_adj_rlx)
!------------------------------------------------------------------------------
! Relax temperature and humidity increments according to the conv. time scale.
!
! This routine modifies the amount by which potential temperature and
! specific humidity are increased to account for convective mixing
! timescale and model timestep.
!
! In cases where moist convection has happened, this routine is
! called twice; once for theta_adjust, and once for theta_noq_adjust.
!------------------------------------------------------------------------------

IMPLICIT NONE

! Inputs
REAL, INTENT(IN) :: theta_start(:)
REAL, INTENT(IN) :: theta_adjust(:)
REAL, INTENT(IN) :: q_start(:)
REAL, INTENT(IN) :: q_adjust(:)
REAL, INTENT(IN) :: timestep

! Outputs
REAL, INTENT(IN OUT) :: theta_adj_rlx(:)
REAL, INTENT(IN OUT) :: q_adj_rlx(:)

! Locals
REAL :: exp_fac

! Can be simplified to
! A_rlx = A_adj - exp() * (A_adj - A_start)

exp_fac = 1.0 - EXP(-timestep / llcs_timescale)

theta_adj_rlx = theta_start + (theta_adjust - theta_start) * exp_fac
q_adj_rlx = q_start + (q_adjust - q_start) * exp_fac

END SUBROUTINE time_relax


SUBROUTINE conserve_moisture(                                                  &
    ! Inputs
    nlevels, p_full, delta_p, q_start, q_adj_rlx,                              &
    ! Outputs
    q_conserve)

!------------------------------------------------------------------------------
! Conserve moisture after convection adjustment.
!
! This routine ensures that at this point in the scheme there is no net
! change in moisture within the column (this will change later, as some
! is lost through rainfall).
!
! To do this, the routine find_bounds is used to find the top and
! bottom of convective events that have occured within the column.
!
! For each event, we find the total mass of water contained on the levels
! involved in this event both before and after convection. We then take
! the ratio of water : before/after and multiply the specfic humidity on
! each level involved in the event by this ratio.
!------------------------------------------------------------------------------

IMPLICIT NONE

! Inputs
INTEGER, INTENT(IN) :: nlevels

REAL, INTENT(IN) :: p_full   (1:nlevels)
REAL, INTENT(IN) :: delta_p  (1:nlevels)
REAL, INTENT(IN) :: q_start  (1:nlevels)
REAL, INTENT(IN) :: q_adj_rlx(1:nlevels)

! Outputs
REAL, INTENT(OUT) :: q_conserve(1:nlevels)

! Locals
INTEGER :: bottom ! used by find_bounds()
INTEGER :: top ! used by find_bounds()
INTEGER :: j_lev
INTEGER :: k
INTEGER :: n

REAL :: q_start_total_mass
REAL :: q_adj_rlx_total_mass
REAL :: q_inc_ratio
REAL :: q_start_mass  (1:nlevels)
REAL :: q_adj_rlx_mass(1:nlevels)
REAL :: q_inc         (1:nlevels)

LOGICAL :: stopnow ! used by find_bounds()

stopnow = .FALSE.
q_conserve = q_adj_rlx

! Find increase in specific humidity after convection
! find_bounds() uses the difference between the start and adjusted profile
! to determine where a convective event begins and ends. For logistical
! reasons it requires a positive difference. For more information, see
! description of find_bounds().
q_inc = ABS(q_adj_rlx - q_start)

n = 2 ! part of find_bounds set up.

DO j_lev = 1, nlevels
  q_start_mass(:) = 0.0
  q_adj_rlx_mass(:) = 0.0

  IF (j_lev == n) THEN
    CALL find_bounds(                                                          &
        ! Inputs
        nlevels, q_inc,                                                        &
        ! Inputs/outputs
        n, stopnow,                                                            &
        ! Outputs
        bottom, top)

    IF (.NOT. stopnow) THEN
      DO k = bottom, top ! find masses on each level within event
        q_start_mass(k) = ABS(q_start(k)) * delta_p(k) / grav
        q_adj_rlx_mass(k) = ABS(q_adj_rlx(k)) * delta_p(k) / grav
      END DO

      ! find total masses for each event
      q_start_total_mass = SUM(q_start_mass(bottom:top))
      q_adj_rlx_total_mass = SUM(q_adj_rlx_mass(bottom:top))

      ! find ratio before / after
      q_inc_ratio = q_start_total_mass / q_adj_rlx_total_mass

      ! set q_conserve accordingly
      DO k = bottom, top
        q_conserve(k) = ABS(q_adj_rlx(k)) * q_inc_ratio
      END DO

    ELSE
      EXIT
    END IF
  END IF

END DO

END SUBROUTINE conserve_moisture


SUBROUTINE rain(                                                               &
    ! Inputs
    nlevels, theta_start, theta_noq_adj_rlx, q_conserve, qcl_start, xmin,      &
    p_full, delta_p, timestep,                                                  &
    ! Inputs/outputs
    theta_adj_rlx,                                                             &
    ! Outputs
    rain_sec, q_final, qcl_inc, cloud_frac)
!------------------------------------------------------------------------------
! Calculate latent heating and rain out the equivalent amount of moisture.
!
! This routine calculates the total amount of latent heating within a
! convective event that was required for the moist convection to take
! place.
!
! It then removes the mass of water necessary to do this in the form of
! instantaneous precipitation and adjusts the specific humidity of the
! levels under consideration accordingly.
!
! Total rainfall in the column is found via Q = mass * L where Q is
! latent heating.
!------------------------------------------------------------------------------

IMPLICIT NONE

! Inputs
INTEGER, INTENT(IN) :: nlevels

REAL, INTENT(IN) :: theta_start      (1:nlevels)
REAL, INTENT(IN) :: theta_noq_adj_rlx(1:nlevels)
REAL, INTENT(IN) :: q_conserve       (1:nlevels)
REAL, INTENT(IN) :: qcl_start        (1:nlevels)
REAL, INTENT(IN) :: xmin             (1:nlevels)
REAL, INTENT(IN) :: p_full           (1:nlevels)
REAL, INTENT(IN) :: delta_p          (1:nlevels)
REAL, INTENT(IN) :: timestep

! Inputs/Outputs
REAL, INTENT(IN OUT) :: theta_adj_rlx (1:nlevels)

! Outputs
REAL, INTENT(OUT) :: rain_sec
REAL, INTENT(OUT) :: q_final         (1:nlevels)
REAL, INTENT(OUT) :: qcl_inc         (1:nlevels)
REAL, INTENT(OUT) :: cloud_frac      (1:nlevels)

! Locals
INTEGER :: j_lev
INTEGER :: k
INTEGER :: n
INTEGER :: bottom ! used by find_bounds()
INTEGER :: top ! used by find_bounds()

REAL :: event_q_total
REAL :: total_latent
REAL :: total_latent_new
REAL :: rain_out
REAL :: rain_out_new
REAL :: latent_ratio
REAL :: temp_new
REAL :: q_level(1:nlevels)
REAL :: latent(1:nlevels)
REAL :: temp(1:nlevels)
REAL :: temp_noq(1:nlevels)
REAL :: theta_inc(1:nlevels)
REAL :: qcl_excess(1:nlevels)
REAL :: qcl_tot

LOGICAL :: stopnow ! used by find_bounds()
LOGICAL :: is_conv(1:nlevels)

stopnow = .FALSE. ! initialise stop now
is_conv(:) = .FALSE.
rain_sec = 0.0

q_final = q_conserve

! Find increase in potential temp after convection
! find_bounds uses the difference between the start and adjusted profile
! (as before, see description in moisture_conserve)
theta_inc = ABS(theta_adj_rlx - theta_start)
qcl_inc(:) = 0.0
qcl_excess(:) = 0.0
cloud_frac(:) = 0.0
qcl_tot = 0.0

n = 2 ! part of find_bounds set up.

! calculate temperature for dummy 'only dry convection' theta profile
CALL calc_temp(theta_noq_adj_rlx, p_full, temp_noq)

DO j_lev = 1, nlevels
  IF (j_lev == n) THEN
    CALL find_bounds(                                                          &
        ! Inputs
        nlevels, theta_inc,                                                    &
        ! Inputs/outputs
        n, stopnow,                                                            &
        ! Outputs
        bottom, top)

    IF (stopnow) THEN
      EXIT
    ELSE
      ! calculate temperature array for actual theta profile
      CALL calc_temp(theta_adj_rlx, p_full, temp)

      ! reset all arrays and totals to 0 for new event
      latent(:) = 0.0
      q_level(:) = 0.0
      event_q_total = 0.0
      total_latent = 0.0
      total_latent_new = 0.0
      latent_ratio = 0.0
      rain_out = 0.0
      rain_out_new = 0.0

      ! Rainfall / latent heating adjustment
      DO k = bottom, top
        latent(k) = delta_p(k) / grav * cp_air * (temp(k) - temp_noq(k))
        is_conv(k) = .TRUE.
      END DO

      ! Calculate rainfall within this column from the amount
      ! of latent heating that has taken place.
      total_latent = SUM(latent(bottom:top))
      IF (total_latent > 0.0) THEN
        DO k = bottom, top
          qcl_inc(k) = latent(k) * (grav / delta_p(k)) / hlv
          ! CHANGE HERE WAS HLV*TIMESTEP
          IF (qcl_inc(k) > 0.0) THEN
            ! Cloud fraction representing detrainment
            ! Needs to be high to emulate the effect of anvils. At the same
            ! time, because the current scheme implies that rain falls "below"
            ! cloud fraction, large values of cloud_frac lead to
            ! precipitation cores being spread out across the area of anvils,
            ! thus making precipitation too slow and evaporating too much.
            cloud_frac(k) = llcs_detrain_coef * timestep / llcs_timescale
          END IF
        END DO

        rain_out = total_latent / hlv

        ! find total water on levels involved in convective event.
        DO k = bottom, top
          q_level(k) = (delta_p(k) / grav) * q_conserve(k)
        END DO
        event_q_total = SUM(q_level(bottom:top))

        ! The following code handles situations where the amount of rain
        ! demanded by a convective event is greater than the amount of
        ! rain stored in the levels participating in the convective event.
        ! In this scenario we limit rain_out to the amount of moisture
        ! present in the event and revise our convective potential
        ! temperature adjustment accordingly.
        IF (rain_out > event_q_total) THEN
          rain_out_new = event_q_total
          total_latent_new = rain_out_new * hlv
          latent_ratio  = total_latent_new / total_latent
          latent = latent * latent_ratio
          DO k = bottom, top
            temp_new = temp_noq(k) + latent(k) / ((delta_p(k) / grav) * cp_air)
            theta_adj_rlx(k) = temp_new * (p_zero/p_full(k))**kappa_use 
          END DO
          qcl_inc = qcl_inc * event_q_total / rain_out
          event_q_total = rain_out_new
          rain_out = rain_out_new
        END IF

        ! Rescale specific humidity on levels within this event to
        ! account for rainfall
        DO k = bottom, top
          q_level(k) = q_level(k) * (event_q_total - rain_out) / event_q_total
          q_final(k) = (grav/delta_p(k)) * q_level(k)
        END DO
      END IF ! total_latent > 0
    END IF ! stopnow
  END IF ! j_lev == n
END DO ! j_lev-loop

SELECT CASE (llcs_cloud_precip)
CASE (llcs_opt_all_rain)
  ! produce rainfall
  rain_sec = SUM(qcl_inc * delta_p / grav) / timestep
  ! zero the cloud increments
  qcl_inc(:) = 0.0
  cloud_frac(:) = 0.0

CASE (llcs_opt_all_cloud)
  ! zero the rain-rate and hand all condensate to LS cloud
  rain_sec = 0.0

CASE (llcs_opt_crit_condens)
  ! Rain out excess of moisture and pass the rest of condensate to LS cloud
  DO k = 1, nlevels
    qcl_tot = qcl_start(k) + qcl_inc(k)
    IF ( (qcl_tot > xmin(k)) .AND. is_conv(k) ) THEN
      ! rain_out = (qcl_tot - xmin(k)) * delta_p(k) / g
      ! Excess of qcl to be rained out
      qcl_excess(k) = qcl_tot - xmin(k)
      ! The next line is equivalent to setting qcl_tot = xmin
      qcl_inc(k) = xmin(k) - qcl_start(k)
    END IF
  END DO
  rain_sec = SUM(qcl_excess * delta_p / grav) / timestep

CASE (llcs_opt_const_frac)
  ! Fixed fraction of rain vs cloud
  rain_sec = SUM(qcl_inc * llcs_rain_frac * delta_p / grav) / timestep
  qcl_inc(:) = qcl_inc(:) * (1 - llcs_rain_frac)
  cloud_frac(:) = cloud_frac(:) * (1 - llcs_rain_frac)
END SELECT

END SUBROUTINE rain


SUBROUTINE conserve_enthalpy(                                                  &
    ! Inputs
    nlevels, theta_start, theta_adj_rlx, theta_noq_adj_rlx, p_full, delta_p,   &
    ! Outputs
    temp_final)

!------------------------------------------------------------------------------
! Conserve energy that was added to the column in the form of dry heating.
!
! Within a given convective event between the levels bottom and top
! (found by find_bounds), we find out how much energy was required to
! allow the increase in potential temperature from theta_start to
! theta_noq_adj_rlx (our dummy dry profile).
!
! In order to calculate how much energy we have added during the event
! under consideration we use dQ = cp_air * dT * mass.
!------------------------------------------------------------------------------

IMPLICIT NONE

! Inputs
INTEGER, INTENT(IN) :: nlevels

REAL, INTENT(IN) :: theta_start      (1:nlevels)
REAL, INTENT(IN) :: theta_adj_rlx    (1:nlevels)
REAL, INTENT(IN) :: theta_noq_adj_rlx(1:nlevels)
REAL, INTENT(IN) :: p_full           (1:nlevels)
REAL, INTENT(IN) :: delta_p          (1:nlevels)

! Outputs
REAL, INTENT(OUT) :: temp_final      (1:nlevels)

! Locals
INTEGER :: bottom ! used by find_bounds()
INTEGER :: top ! used by find_bounds()
INTEGER :: n
INTEGER :: j_lev
INTEGER :: k

REAL :: total_heating
REAL :: total_mass
REAL :: temp_correction
REAL :: temp_start      (1:nlevels)
REAL :: temp_adj_rlx    (1:nlevels)
REAL :: temp_noq_adj_rlx(1:nlevels)
REAL :: dry_heating     (1:nlevels)
REAL :: mass            (1:nlevels)
REAL :: theta_dry_inc   (1:nlevels)

LOGICAL :: stopnow ! used by find_bounds()

! Find temperatures
CALL calc_temp(theta_start, p_full, temp_start)
CALL calc_temp(theta_adj_rlx, p_full, temp_adj_rlx)
CALL calc_temp(theta_noq_adj_rlx, p_full, temp_noq_adj_rlx)
temp_final = temp_adj_rlx

stopnow = .FALSE.

theta_dry_inc = theta_noq_adj_rlx - theta_start

n = 2 ! required for find_bounds

DO j_lev = 1, nlevels
  IF (j_lev == n) THEN
    CALL find_bounds(                                                          &
        ! Inputs
        nlevels, theta_dry_inc,                                                &
        ! Inputs/outputs
        n, stopnow,                                                            &
        ! Outputs
        bottom, top)
    IF (.NOT. stopnow) THEN
      DO k = bottom, top
        mass(k) = 0.0
        dry_heating(k) = 0.0
        mass(k) = delta_p(k) / grav
        dry_heating(k) = cp_air * (temp_noq_adj_rlx(k) - temp_start(k)) * mass(k)
      END DO
      total_mass = SUM(mass(bottom:top))
      total_heating = SUM(dry_heating(bottom:top))
      temp_correction = total_heating / (cp_air * total_mass)
      DO k = bottom, top
        temp_final(k) = temp_adj_rlx(k) - temp_correction
      END DO
    ELSE
      EXIT
    END IF
  END IF
END DO

END SUBROUTINE conserve_enthalpy


SUBROUTINE critical_condensate(                                                &
      ! Inputs
      nlevels, theta_start, q_start, p_full,                                   &
      ! Inputs/outputs
      xmin)
!------------------------------------------------------------------------------
! Calculate critical cloud condensate from initial temperature and humidity,
! given the free parameters: minimum, maximum and qsat scaling.
!
! The result is used as a threshold to separate convective rainfall and cloud
! condensed water (the latter is then passed to the cloud scheme)
! The subroutine is called if llcs_cloud_precip is equal to 
! llcs_opt_crit_condens.
!------------------------------------------------------------------------------


IMPLICIT NONE

! Inputs
INTEGER, INTENT(IN) :: nlevels

REAL, INTENT(IN) :: theta_start(1:nlevels)
REAL, INTENT(IN) :: q_start    (1:nlevels)
REAL, INTENT(IN) :: p_full     (1:nlevels)

! Inputs/outputs
REAL, INTENT(IN OUT) :: xmin    (1:nlevels)

! Locals
INTEGER :: j_lev

REAL :: temp
REAL :: q_sat

DO j_lev = 1, nlevels
  temp = theta_start(j_lev) * (p_full(j_lev)/p_zero) ** (kappa_use)
  CALL qsat_calc(p_full(j_lev), temp, q_sat)

  ! Equivalent to opt. 4 in cloud_w_6a scheme in Met Office UM 
  xmin(j_lev) = MAX(MIN(mparwtr, fac_qsat * q_sat), qlmin)
END DO

END SUBROUTINE critical_condensate

! ----------------------------------------------------------------------
! -- END OF ROUTINES CALLED BY CONVECTION SCHEME -- !
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! -- EXTRA ROUTINES -- !
! ----------------------------------------------------------------------

SUBROUTINE find_bounds(                                                        &
    ! Inputs
    nlevels, inc,                                                              &
    ! Inputs/outputs
    n, stopnow,                                                                &
    ! Outputs
    bottom, top)
!------------------------------------------------------------------------------
! FIND BOUNDS
! This routine finds the bounds for a convective event from either the
! potential temperature or specific humidity increment for that event.
!
! Note that this routine also outputs the integer n which is the level
! above the top level for the convective event under consideration. In
! other words, this is the first level above a convective event where a
! different convective event could be found.
!------------------------------------------------------------------------------

IMPLICIT NONE

! Inputs
INTEGER, INTENT(IN) :: nlevels

REAL, INTENT(IN) :: inc(1:nlevels)

! Inputs/Outputs
INTEGER, INTENT(IN OUT) :: n

LOGICAL, INTENT(IN OUT) :: stopnow

! Outputs
INTEGER, INTENT(OUT) :: bottom, top

! Locals
INTEGER :: j_lev

bottom = -1

DO j_lev = n, nlevels
  ! find bottom
  IF (bottom == -1) THEN
    IF (inc(j_lev) > 0) THEN
      bottom = j_lev - 1
    ELSE IF (j_lev == nlevels) THEN
      stopnow = .TRUE.
      EXIT
    END IF
  END IF

  ! find top
  IF (j_lev == nlevels) THEN
    top = j_lev - 1
    EXIT
  ELSE IF ( (j_lev - 1 /= bottom) .AND. (inc(j_lev) == 0) ) THEN
    IF (bottom /= -1) THEN
      top = j_lev - 1
      EXIT
    END IF
  END IF

END DO

n = top + 1

END SUBROUTINE find_bounds

! ----------------------------------------------------------------------

SUBROUTINE calc_temp(theta, p_full, temp)

! CALCULATE TEMPERATURE
! This routine calculates and array containing temperatures from an
! array containing potential temperatures.

IMPLICIT NONE

REAL, INTENT(IN) :: theta(:), p_full(:)
REAL, INTENT(OUT) :: temp(:)

temp = theta * (p_full/p_zero)**(kappa_use)

END SUBROUTINE calc_temp

! ----------------------------------------------------------------------

SUBROUTINE calc_theta(temp, p_full, theta)

! CALCULATE TEMPERATURE
! This routine calculates and array containing temperatures from an
! array containing potential temperatures.

IMPLICIT NONE

REAL, INTENT(IN) :: temp(:), p_full(:)
REAL, INTENT(OUT) :: theta(:)

theta = temp * (p_zero/p_full)**(kappa_use)

END SUBROUTINE calc_theta

subroutine qsat_calc(pfull,T,qsat) !NTL subroutine modified from 2006 FMS MoistModel file moist_processes.f90 (v14 2012/06/22 14:50:00).

        IMPLICIT NONE


        REAL, INTENT (IN)    :: pfull, T
        REAL, INTENT (OUT)   :: qsat

        REAL :: esat

!-----------------------------------------------------------------------
!       Calculate saturation specific humidity 
!       This is calculated according to the formula:
!
!       qsat   = (epsilon*esat) / [pfull  -  (1.-epsilon)*esat] 
!
!       Where epsilon = RDGAS/RVGAS is 'repsilon' in the code 
!
!       and where 1- epsilon is 'one_minus_epsilon' in the code 
!
!       Note that qsat does not have its proper value
!       until all of the following code has been executed.  That
!       is, qsat is used to store intermediary results
!       in forming the full solution.

        !calculate water saturated vapor pressure from table
        !and store temporarily in the variable esat
        CALL escomp(T,esat)						!same as escomp

        !calculate denominator in qsat formula
        if(do_simple) then
          qsat = pfull
        else
          qsat = pfull-one_minus_epsilon*esat
        endif

        !limit denominator to esat, and thus qs to epsilon
        !this is done to avoid blow up in the upper stratosphere
        !where pfull ~ esat
        qsat= MAX(qsat,esat)

        !calculate qsat
        qsat = repsilon * esat / qsat
END SUBROUTINE qsat_calc

subroutine rsat_calc(pfull,T,rsat) 

        IMPLICIT NONE


        REAL, INTENT (IN)  :: pfull, T
        REAL, INTENT (OUT) :: rsat

        REAL :: esat

!-----------------------------------------------------------------------
!       Calculate saturation mixing ratio  
!       This is calculated according to the formula:
!
!       rsat   = (epsilon*esat) / [pfull - esat] 
!
!       Where epsilon = RDGAS/RVGAS is 'repsilon' in the code 
!
!       and where 1- epsilon is 'one_minus_epsilon' in the code 
!
!       Note that rsat does not have its proper value
!       until all of the following code has been executed.  That
!       is, rsat is used to store intermediary results
!       in forming the full solution.

        !calculate water saturated vapor pressure from table
        !and store temporarily in the variable esat
        CALL escomp(T,esat)						!same as escomp

        !calculate denominator in rsat formula
        rsat = pfull - esat
        ! WARNING: rsat may get large as pfull ~ esat in stratosphere. Unlike for qsat 
        ! calculation, there is no 'do simple' fix for this. 

        !calculate qsat
        rsat = repsilon * esat / rsat
END SUBROUTINE rsat_calc



! ----------------------------------------------------------------------
! -- END OF EXTRA ROUTINES -- !
! ----------------------------------------------------------------------


END MODULE llcs_mod
