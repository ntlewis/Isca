

module titan_gray_rad_mod

! ==================================================================================
! ==================================================================================

   use fms_mod,               only: open_file, check_nml_error, &
                                    mpp_pe, close_file, error_mesg, &
                                    NOTE, FATAL,  uppercase

   use constants_mod,         only: stefan, cp_air, grav, pstd_mks, pstd_mks_earth, seconds_per_sol, orbital_period

   use    diag_manager_mod,   only: register_diag_field, send_data

   use    time_manager_mod,   only: time_type, length_of_year, length_of_day, &
                                    operator(+), operator(-), operator(/=), get_time, get_date

   use astronomy_mod,         only: astronomy_init, diurnal_solar

   use interpolator_mod,      only: interpolate_type, interpolator_init, interpolator, ZERO, interpolator_end

!==================================================================================
implicit none
private
!==================================================================================

! version information

character(len=128) :: version='$Id: titan_gray_rad.F90,v 1.1.2.1 2021/02/18$'
character(len=128) :: tag='Two-stream gray atmosphere for titan'

! Version Details
! [2021/02/18] <Neil Lewis>:  Began coding 

! public interfaces

public :: titan_gray_rad_init, titan_gray_rad_down, titan_gray_rad_up, titan_gray_rad_end
!==================================================================================


! module variables
logical :: initialized     = .false.

! varcp params
logical :: var_cp = .false. 


! optical depth parameters ! defaults to McKay 
real    :: tau_zero = 4.22 ! something 
real    :: pres_zero = 150000 ! something 
real    :: linear_tau = 0.0
real    :: nonlin_tau_exponent = 4./3.
real    :: lin_tau_exponent = 1.
real    :: sw_tau_exponent = 1.
logical :: do_rc_optical_depth = .false. 
real    :: tau_rc = -1  ! only used if do_rc_optical_depth = .true. 
real    :: pres_rc = -1 ! should be set by user

! sw parameters ! defaults to McKay 
real    :: bond_albedo = 0.27
real    :: sw_k1 = 35.167 ! fix these
real    :: sw_k2 = 0.0613 
real    :: sw_k3 = 0.0
real    :: sw_gamma = 0.0
real    :: sw_k1_frac = 4./7.
logical :: do_sw_split = .true. 

logical :: do_sw_venus = .false.
logical :: do_log_sw = .false. 
logical :: do_log_sw2 = .false.
logical :: do_log_sw3 = .false. 
real    :: logtau_B = 0.09229!0.272 / 2.646
real    :: logtau_sw_s  = 1.897
real    :: logtau_0 = 356.7
real    :: logtau_siglev = 0.0002866
real    :: logtau_pres_zero = 150000
real    :: logf = 1./3.

! insolation parameters 
real    :: solar_constant  = 3.6*4.
real    :: del_sol         = 1.4
real    :: del_sw          = 0.0
logical :: do_seasonal     = .false.
integer :: solday          = -10  !s Day of year to run perpetually if do_seasonal=True and solday>0
real    :: equinox_day     = 0.75 !s Fraction of year [0,1] where NH autumn equinox occurs (only really useful if calendar has defined months).
logical :: use_time_average_coszen = .false. !s if .true., then time-averaging is done on coszen so that insolation doesn't depend on timestep
real    :: dt_rad_avg     = -1
logical :: simple_diurnal = .false.
logical :: tidally_locked = .false. 
logical :: perp_eq = .false.
real    :: noon_longitude_tl = 180.
real    :: diurnal_period = -1. !diurnal period in earth days. -ve yields east to west 



real :: diabatic_acce = 1.


real, allocatable, dimension(:,:)   :: insolation, p2, b_surf
real, allocatable, dimension(:,:,:) :: b, tdt_rad, tdt_solar
real, allocatable, dimension(:,:,:) :: lw_up, lw_down, lw_flux, sw_up, sw_down, sw_flux, rad_flux
real, allocatable, dimension(:,:,:) :: tau, del_tau 
real, allocatable, dimension(:,:)   :: olr, net_lw_surf, toa_sw_in, coszen, fracsun

real, save :: pi, deg_to_rad , rad_to_deg




namelist/titan_gray_rad_nml/ tau_rc, tau_zero, pres_rc, pres_zero, & 
                             do_rc_optical_depth, linear_tau, &
                             nonlin_tau_exponent, &
                             bond_albedo, sw_k1, sw_k2, sw_k1_frac, & 
                             sw_k3, sw_gamma, &
                             do_sw_split, solar_constant, del_sol, del_sw, &
                             do_seasonal, solday, equinox_day, dt_rad_avg, & 
                             use_time_average_coszen, diabatic_acce, & 
                             noon_longitude_tl, tidally_locked, &
                             simple_diurnal, diurnal_period, perp_eq, &
                             do_log_sw, logtau_sw_s, logtau_B, logtau_0, &
                             logtau_siglev,logtau_pres_zero, do_log_sw2, logf, &
                             do_log_sw3, lin_tau_exponent, sw_tau_exponent, var_cp, do_sw_venus

!==================================================================================
!-------------------- diagnostics fields -------------------------------

integer :: id_olr, id_swdn_sfc, id_swdn_toa, id_net_lw_surf, id_lwdn_sfc, id_lwup_sfc, &
           id_tdt_rad, id_tdt_solar, id_flux_rad, id_flux_lw, id_flux_sw, id_coszen, id_fracsun, id_tau 

character(len=14), parameter :: mod_name = 'titan_gray_rad'

real :: missing_value = -999.


contains


! ==================================================================================
! ==================================================================================


subroutine titan_gray_rad_init(is, ie, js, je, num_levels, axes, Time, lonb, latb, dt_real)

!-------------------------------------------------------------------------------------
integer, intent(in), dimension(4) :: axes
type(time_type), intent(in)       :: Time
integer, intent(in)               :: is, ie, js, je, num_levels
real ,dimension(:,:),intent(in),optional :: lonb,latb !s Changed to 2d arrays as 2013 interpolator expects this.
real, intent(in)                  :: dt_real !s atmospheric timestep, used for radiation averaging
!-------------------------------------------------------------------------------------
integer, dimension(3) :: half = (/1,2,4/)
integer :: ierr, io, unit
!-----------------------------------------------------------------------------------------
! read namelist and copy to logfile

unit = open_file ('input.nml', action='read')
ierr=1
do while (ierr /= 0)
   read  (unit, nml=titan_gray_rad_nml, iostat=io, end=10)
   ierr = check_nml_error (io, 'titan_gray_rad_nml')
enddo
10 call close_file (unit)

unit = open_file ('logfile.out', action='append')
if ( mpp_pe() == 0 ) then
  write (unit,'(/,80("="),/(a))') trim(version), trim(tag)
  write (unit, nml=titan_gray_rad_nml)
endif
call close_file (unit)

pi         = 4. * atan(1.)
deg_to_rad = pi/180.
rad_to_deg = 180./pi

call astronomy_init

if(dt_rad_avg .le. 0) dt_rad_avg = dt_real !s if dt_rad_avg is set to a value in nml then it will be used instead of dt_real









initialized = .true.

allocate (b                (ie-is+1, je-js+1, num_levels))
allocate (tdt_rad          (ie-is+1, je-js+1, num_levels))
allocate (tdt_solar        (ie-is+1, je-js+1, num_levels))


allocate (tau              (ie-is+1, je-js+1, num_levels+1))
allocate (del_tau          (ie-is+1, je-js+1, num_levels))
allocate (lw_up            (ie-is+1, je-js+1, num_levels+1))
allocate (lw_down          (ie-is+1, je-js+1, num_levels+1))
allocate (lw_flux          (ie-is+1, je-js+1, num_levels+1))
allocate (sw_up            (ie-is+1, je-js+1, num_levels+1))
allocate (sw_down          (ie-is+1, je-js+1, num_levels+1))
allocate (sw_flux          (ie-is+1, je-js+1, num_levels+1))
allocate (rad_flux         (ie-is+1, je-js+1, num_levels+1))

allocate (b_surf           (ie-is+1, je-js+1))
allocate (olr              (ie-is+1, je-js+1))
allocate (net_lw_surf      (ie-is+1, je-js+1))
allocate (toa_sw_in        (ie-is+1, je-js+1))

allocate (insolation       (ie-is+1, je-js+1))
allocate (p2               (ie-is+1, je-js+1))
allocate (coszen           (ie-is+1, je-js+1))
allocate (fracsun          (ie-is+1, je-js+1)) !jp from astronomy.f90 : fraction of sun on surface



!-----------------------------------------------------------------------
!------------ initialize diagnostic fields ---------------

    id_olr = &
    register_diag_field ( mod_name, 'olr', axes(1:2), Time, &
               'outgoing longwave radiation', &
               'W/m2', missing_value=missing_value               )
    id_swdn_sfc = &
    register_diag_field ( mod_name, 'swdn_sfc', axes(1:2), Time, &
               'Absorbed SW at surface', &
               'W/m2', missing_value=missing_value               )
    id_swdn_toa = &
    register_diag_field ( mod_name, 'swdn_toa', axes(1:2), Time, &
               'SW flux down at TOA', &
               'W/m2', missing_value=missing_value               )
    id_lwup_sfc = &
    register_diag_field ( mod_name, 'lwup_sfc', axes(1:2), Time, &
               'LW flux up at surface', &
               'W/m2', missing_value=missing_value               )

    id_lwdn_sfc = &
    register_diag_field ( mod_name, 'lwdn_sfc', axes(1:2), Time, &
               'LW flux down at surface', &
               'W/m2', missing_value=missing_value               )

    id_net_lw_surf = &
    register_diag_field ( mod_name, 'net_lw_surf', axes(1:2), Time, &
               'Net upward LW flux at surface', &
               'W/m2', missing_value=missing_value               )

    id_tdt_rad = &
        register_diag_field ( mod_name, 'tdt_rad', axes(1:3), Time, &
               'Temperature tendency due to radiation', &
               'K/s', missing_value=missing_value               )

    id_tdt_solar = &
        register_diag_field ( mod_name, 'tdt_solar', axes(1:3), Time, &
               'Temperature tendency due to solar radiation', &
               'K/s', missing_value=missing_value               )

    id_flux_rad = &
        register_diag_field ( mod_name, 'flux_rad', axes(half), Time, &
               'Total radiative flux (positive up)', &
               'W/m^2', missing_value=missing_value               )
    id_flux_lw = &
        register_diag_field ( mod_name, 'flux_lw', axes(half), Time, &
               'Net longwave radiative flux (positive up)', &
               'W/m^2', missing_value=missing_value               )
    id_flux_sw = &
        register_diag_field ( mod_name, 'flux_sw', axes(half), Time, &
               'Net shortwave radiative flux (positive up)', &
               'W/m^2', missing_value=missing_value               )

    id_coszen  = &
               register_diag_field ( mod_name, 'coszen', axes(1:2), Time, &
                 'cosine of zenith angle', &
                 'none', missing_value=missing_value      )
    id_fracsun  = &
               register_diag_field ( mod_name, 'fracsun', axes(1:2), Time, &
                 'daylight fraction of time interval', &
                 'none', missing_value=missing_value      )

    id_tau      = &
                 register_diag_field ( mod_name, 'tau', axes(half), Time, &
                   'tau', &
                   'none', missing_value=missing_value      )


return
end subroutine titan_gray_rad_init

! ==================================================================================

subroutine titan_gray_rad_down (is, js, Time_diag, lat, lon, p_full, p_half, t,         &
                           net_surf_sw_down, surf_lw_down, albedo, q)

! Begin the radiation calculation by computing downward fluxes.
! This part of the calculation does not depend on the surface temperature.

integer, intent(in)                 :: is, js
type(time_type), intent(in)         :: Time_diag
real, intent(in), dimension(:,:)    :: lat, lon, albedo
real, intent(out), dimension(:,:)   :: net_surf_sw_down
real, intent(out), dimension(:,:)   :: surf_lw_down
real, intent(in), dimension(:,:,:)  :: t, q,  p_half, p_full
integer :: i, j, k, n, dyofyr

integer :: seconds, year_in_s, days, months, years, hours, minutes 
real :: r_seconds, frac_of_day, frac_of_year, gmt, time_since_ae, rrsun, day_in_s, r_solday, r_total_seconds, r_days, r_dt_rad_avg, dt_rad_radians, noon_longitude, r_months, r_years, r_hours, r_minutes
logical :: used

!logical :: no_exp 
real, dimension(size(t,1),size(t,2)) :: exp_insol

real ,dimension(size(q,1),size(q,2),size(q,3)) :: co2f
real, dimension(size(t,1),size(t,2),size(t,3)+1) :: sw_down_k1, sw_down_k2, sw_down_k3



n = size(t,3)


! albedo(:,:) = albedo_value !s albedo now set in mixed_layer_init.

! zero fluxes before starting 
sw_down_k1 = 0.0 
sw_down_k2 = 0.0 
sw_down_k3 = 0.0 
sw_down    = 0.0 
lw_down    = 0.0 

! =================================================================================
! SHORTWAVE RADIATION

! insolation at TOA
if (do_seasonal) then
  ! Seasonal Cycle: Use astronomical parameters to calculate insolation
  call get_time(Time_diag, seconds, days)
  call get_time(length_of_year(), year_in_s)
  r_seconds = real(seconds)
  day_in_s = length_of_day()
  !write(*,*) 'day:', day_in_s 
  !write(*,*) 'year:', year_in_s 
  frac_of_day = r_seconds / day_in_s

  if(solday .ge. 0) then
      r_solday=real(solday)
      frac_of_year = (r_solday*day_in_s) / year_in_s
  else
      r_days=real(days)
      r_total_seconds=r_seconds+(r_days*86400.)
      frac_of_year = r_total_seconds / year_in_s
  endif
  !write(*,*) 'days:', r_days 
  !write(*,*) 'frac:', frac_of_year

  gmt = abs(mod(frac_of_day, 1.0)) * 2.0 * pi

  time_since_ae = modulo(frac_of_year-equinox_day, 1.0) * 2.0 * pi

  if(use_time_average_coszen) then
     r_dt_rad_avg=real(dt_rad_avg)
     dt_rad_radians = 2.0 * pi !(r_dt_rad_avg/day_in_s)*2.0*pi
     call diurnal_solar(lat, lon, gmt, time_since_ae, coszen, fracsun, rrsun, dt_rad_radians)
  else
     call diurnal_solar(lat, lon, gmt, time_since_ae, coszen, fracsun, rrsun)
  end if

elseif (simple_diurnal) then 
   if (tidally_locked) then 
      noon_longitude = noon_longitude_tl*deg_to_rad
   else 
      call get_date(Time_diag, years, months, days, hours, minutes, seconds)
      r_seconds = real(seconds)
      r_minutes = real(minutes)
      r_hours = real(hours)
      r_days = real(days)
      r_months = real(months)
      r_years = real(years)
      r_total_seconds = ((r_years-1.)*360. + (r_months-1.)*30. + (r_days-1.))*86400. + r_hours * 3600. + r_minutes*60. + r_seconds 
      noon_longitude = 2.0 * pi * modulo(r_total_seconds/(diurnal_period*86400.), 1.0)
      !write(*,*) 'years', r_years, 'months', r_months, 'days', r_days, 'hours', r_hours, 'minutes', r_minutes, 'seconds', r_seconds
   endif 
   coszen = cos(lat)*cos(lon - noon_longitude)
   where ((coszen).le.0.0) 
      coszen = 0.0 
   endwhere 

elseif (perp_eq) then 
   coszen = cos(lat) / pi 

else
  ! Default: Averaged Earth insolation at all longitudes
  p2          = (1. - 3.*sin(lat)**2)/4.
  coszen  = 0.25 * (1.0 + del_sol * p2 + del_sw * sin(lat)) 
end if
insolation = solar_constant * (1 - bond_albedo) * coszen 

! compute optical depth 
!del_tau = tau_zero / pres_zero * (linear_tau + nonlin_tau_exponent * (1. - linear_tau) * (p_full(:,:,:) / pres_zero) ** (nonlin_tau_exponent - 1.)) * (p_half(:,:,2:n+1) - p_half(:,:,1:n))

if (do_rc_optical_depth) then 
   where (p_half .lt. pres_rc) 
      tau = tau_rc * (linear_tau * p_half(:,:,:) / pres_rc + (1. - linear_tau) * (p_half(:,:,:) / pres_rc) ** nonlin_tau_exponent)
   elsewhere
      tau = tau_zero * (p_half(:,:,:) / pres_zero) ** nonlin_tau_exponent
   endwhere 
else 
   tau = tau_zero * (linear_tau * (p_half(:,:,:) / pres_zero) ** lin_tau_exponent + (1. - linear_tau) * (p_half(:,:,:) / pres_zero) ** nonlin_tau_exponent)
endif    
tau(:,:,1) = 0.0 


do k = 1, n 
  del_tau(:,:,k) = tau(:,:,k+1) - tau(:,:,k)
enddo 

! compute sw flux 
!no_exp = .false.
if (do_sw_venus) then 
    sw_down(:,:,1) = insolation(:,:)
    !write(*,*) insolation(:,:)
    do k = 1, n
        sw_down(:,:,k+1) =  (60./150.)*insolation(:,:)*exp(-1.5*p_half(:,:,k+1)/p_half(:,:,n+1)) + (1.-(60./150.)) * insolation(:,:)*exp(-500.*p_half(:,:,k+1)/p_half(:,:,n+1))
    enddo 

elseif (do_log_sw) then 
   sw_down(:,:,1) = insolation(:,:)
   do k = 1, n 
      sw_down(:,:,k+1) = min(insolation(:,:), insolation(:,:) * (exp(-logtau_sw_s) - logtau_B * log(p_half(:,:,k+1) / logtau_pres_zero)))   
      exp_insol = insolation(:,:)*exp(-logtau_0*p_half(:,:,k+1) / logtau_pres_zero)
      where (p_half(:,:,k+1)/logtau_pres_zero.lt.logtau_siglev) 
         sw_down(:,:,k+1) = exp_insol 
      endwhere 
      ! if (.not.no_exp) then 
         
      !    if (sw_down(is,js,k+1).gt.exp_insol(is,js)) then 
      !       sw_down(:,:,k+1) = exp_insol 
      !    else 
      !       no_exp=.true. 
      !    endif 
      ! endif 
   enddo 
elseif (do_log_sw2) then 
   sw_down(:,:,1) = insolation(:,:)
   do k = 1, n 
      sw_down(:,:,k+1) = logf*insolation(:,:) * (exp(-logtau_sw_s) - logtau_B * log(p_half(:,:,k+1) / logtau_pres_zero))   
      exp_insol = (1-logf)*insolation(:,:)*exp(-logtau_0*p_half(:,:,k+1) / logtau_pres_zero) + logf*insolation(:,:)
      where (p_half(:,:,k+1)/logtau_pres_zero.lt.logtau_siglev) 
         sw_down(:,:,k+1) = exp_insol 
      endwhere 
   enddo 
elseif (do_log_sw3) then 
      sw_down(:,:,1) = insolation(:,:)
      do k = 1, n 
         sw_down(:,:,k+1) = ((1-logf)*insolation(:,:) * exp(-logtau_0*logtau_siglev**sw_tau_exponent) + logf*insolation(:,:)) * (exp(-logtau_sw_s) - logtau_B * log(p_half(:,:,k+1) / logtau_pres_zero))
         exp_insol = (1-logf)*insolation(:,:)*exp(-logtau_0*(p_half(:,:,k+1) / logtau_pres_zero)**sw_tau_exponent) + logf*insolation(:,:)
         where (p_half(:,:,k+1)/logtau_pres_zero.lt.logtau_siglev) 
            sw_down(:,:,k+1) = exp_insol 
         endwhere 
      enddo 
else    
   if (do_sw_split) then 
   sw_down_k1(:,:,1) = sw_k1_frac * (1. - sw_gamma) * insolation(:,:)
   sw_down_k2(:,:,1) = (1. - sw_k1_frac) * insolation(:,:)
   sw_down_k3(:,:,1) = sw_k1_frac * sw_gamma * insolation(:,:)
   do k = 1, n
      where ((coszen).ne.0.0) 
         sw_down_k1(:,:,k+1) = (2 - sw_k1 / coszen(:,:) * del_tau(:,:,k)) / (2 + sw_k1 / coszen(:,:) * del_tau(:,:,k)) * sw_down_k1(:,:,k)
         sw_down_k2(:,:,k+1) = (2 - sw_k2 / coszen(:,:) * del_tau(:,:,k)) / (2 + sw_k2 / coszen(:,:) * del_tau(:,:,k)) * sw_down_k2(:,:,k)
         sw_down_k3(:,:,k+1) = (2 - sw_k3 / coszen(:,:) * del_tau(:,:,k)) / (2 + sw_k3 / coszen(:,:) * del_tau(:,:,k)) * sw_down_k3(:,:,k)
      endwhere 
   end do
   sw_down = sw_down_k1 + sw_down_k2 + sw_down_k3
   else
      sw_down(:,:,1) = insolation(:,:)
      do k = 1, n 
         where ((coszen).ne.0.0) 
            sw_down(:,:,k+1)   = (2 - sw_k1 / coszen(:,:) * del_tau(:,:,k)) / (2 + sw_k1 / coszen(:,:) * del_tau(:,:,k)) * sw_down(:,:,k)
         endwhere 
      enddo 
   endif 
endif

! compute lw flux
b = stefan*t**4



! compute downward longwave flux by integrating downward
lw_down(:,:,1)      = 0.
do k = 1, n
  lw_down(:,:,k+1) = 2.*b(:,:,k) * del_tau(:,:,k) / (1. + del_tau(:,:,k)) + (1. - del_tau(:,:,k)) / (1. + del_tau(:,:,k)) * lw_down(:,:,k)
end do



! =================================================================================
surf_lw_down     = lw_down(:, :, n+1)
toa_sw_in        = sw_down(:, :, 1)
net_surf_sw_down = sw_down(:, :, n+1) * (1. - albedo)
! =================================================================================




!------- downward lw flux surface -------
if ( id_lwdn_sfc > 0 ) then
   used = send_data ( id_lwdn_sfc, surf_lw_down, Time_diag)
endif
!------- incoming sw flux toa -------
if ( id_swdn_toa > 0 ) then
   used = send_data ( id_swdn_toa, toa_sw_in, Time_diag)
endif
!------- downward sw flux surface -------
if ( id_swdn_sfc > 0 ) then
   used = send_data ( id_swdn_sfc, net_surf_sw_down, Time_diag)
endif

!------- cosine of zenith angle ------------
if ( id_coszen > 0 ) then
   used = send_data ( id_coszen, coszen, Time_diag)
endif

if ( id_tau > 0 ) then
   used = send_data ( id_tau, tau, Time_diag)
endif


return
end subroutine titan_gray_rad_down

! ==================================================================================

subroutine titan_gray_rad_up (is, js, Time_diag, lat, p_half, t_surf, t, tdt, albedo)

! Now complete the radiation calculation by computing the upward and net fluxes.

integer, intent(in)                 :: is, js
type(time_type), intent(in)         :: Time_diag
real, intent(in) , dimension(:,:)   :: lat, albedo
real, intent(in) , dimension(:,:)   :: t_surf
real, intent(in) , dimension(:,:,:) :: t, p_half
real, intent(inout), dimension(:,:,:) :: tdt

real, dimension(size(t,1),size(t,2),size(t,3)) :: cp_air_local

integer :: i, j, k, n

logical :: used


if (var_cp) then 
    cp_air_local(:,:,:) = 1000. * (t(:,:,:) / 460.) ** 0.35 
endif 


n = size(t,3)
b_surf            = stefan*t_surf**4

! zero fluxes before starting 
lw_up = 0.0 
sw_up = 0.0


! compute upward longwave flux by integrating upward

lw_up(:,:,n+1)    = b_surf
do k = n, 1, -1
     lw_up(:,:,k)   = 2. * del_tau(:,:,k) / (1. + del_tau(:,:,k)) * b(:,:,k) + (1. - del_tau(:,:,k)) / (1. + del_tau(:,:,k)) * lw_up(:,:,k+1)
end do


! compute upward shortwave flux (here taken to be constant)
do k = 1, n+1
   where ((coszen).ne.0.0) 
      sw_up(:,:,k)   = albedo(:,:) * sw_down(:,:,n+1)
   endwhere 
end do

! net fluxes (positive up)
lw_flux  = lw_up - lw_down
sw_flux  = sw_up - sw_down
rad_flux = lw_flux + sw_flux

if (.not.var_cp) then 
    do k = 1, n
       tdt_rad(:,:,k)   = diabatic_acce * ( rad_flux(:,:,k+1) - rad_flux(:,:,k) )  &
            * grav/( cp_air*(p_half(:,:,k+1) - p_half(:,:,k)) )

       tdt_solar(:,:,k) = ( sw_flux(:,:,k+1) - sw_flux(:,:,k) )  &
            * grav/( cp_air*(p_half(:,:,k+1) - p_half(:,:,k)) )

       tdt(:,:,k) = tdt(:,:,k) + tdt_rad(:,:,k)
    end do
else 
    do k = 1, n
       tdt_rad(:,:,k)   = diabatic_acce * ( rad_flux(:,:,k+1) - rad_flux(:,:,k) )  &
            * grav/( cp_air_local(:,:,k)*(p_half(:,:,k+1) - p_half(:,:,k)) )

       tdt_solar(:,:,k) = ( sw_flux(:,:,k+1) - sw_flux(:,:,k) )  &
            * grav/( cp_air_local(:,:,k)*(p_half(:,:,k+1) - p_half(:,:,k)) )

       tdt(:,:,k) = tdt(:,:,k) + tdt_rad(:,:,k)
    end do
endif 


olr         = lw_up(:,:,1)
net_lw_surf = lw_flux(:, :, n+1)

!------- outgoing lw flux toa (olr) -------
if ( id_olr > 0 ) then
   used = send_data ( id_olr, olr, Time_diag)
endif
!------- upward lw flux surface -------
if ( id_lwup_sfc > 0 ) then
   used = send_data ( id_lwup_sfc, b_surf, Time_diag)
endif
!------- net upward lw flux surface -------
if ( id_net_lw_surf > 0 ) then
   used = send_data ( id_net_lw_surf, net_lw_surf, Time_diag)
endif
!------- temperature tendency due to radiation ------------
if ( id_tdt_rad > 0 ) then
   used = send_data ( id_tdt_rad, tdt_rad, Time_diag)
endif
!------- temperature tendency due to solar radiation ------------
if ( id_tdt_solar > 0 ) then
   used = send_data ( id_tdt_solar, tdt_solar, Time_diag)
endif
!------- total radiative flux (at half levels) -----------
if ( id_flux_rad > 0 ) then
   used = send_data ( id_flux_rad, rad_flux, Time_diag)
endif
!------- longwave radiative flux (at half levels) --------
if ( id_flux_lw > 0 ) then
   used = send_data ( id_flux_lw, lw_flux, Time_diag)
endif
if ( id_flux_sw > 0 ) then
   used = send_data ( id_flux_sw, sw_flux, Time_diag)
endif


return
end subroutine titan_gray_rad_up

! ==================================================================================


subroutine titan_gray_rad_end

deallocate (b, tdt_rad, tdt_solar)
deallocate (lw_up, lw_down, lw_flux, sw_up, sw_down, sw_flux, rad_flux)
deallocate (b_surf, olr, net_lw_surf, toa_sw_in)
deallocate (tau, del_tau)
deallocate (insolation, p2) !s albedo



end subroutine titan_gray_rad_end

! ==================================================================================

end module titan_gray_rad_mod
