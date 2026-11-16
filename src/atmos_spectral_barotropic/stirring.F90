module stirring_mod

!-----------------------------------------------------------------------
!                   GNU General Public License                        
!                                                                      
! This program is free software; you can redistribute it and/or modify it and  
! are expected to follow the terms of the GNU General Public License  
! as published by the Free Software Foundation; either version 2 of   
! the License, or (at your option) any later version.                 
!                                                                      
! This program is distributed in the hope that it will be useful, but WITHOUT    
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    
! License for more details.                                           
!                                                                      
! For the full text of the GNU General Public License,                
! write to: Free Software Foundation, Inc.,                           
!           675 Mass Ave, Cambridge, MA 02139, USA.                   
! or see:   http://www.gnu.org/licenses/gpl.html                      
!-----------------------------------------------------------------------

! Two stirring options are available, selected via the namelist:
!
! (1) Ornstein-Uhlenbeck (OU) process  [do_white_forcing = .false., default]
!     As described in:
!     Vallis, Gerber, Kushner, Cash, 2004: A Mechanism and Simple Dynamical
!     Model of the North Atlantic Oscillation and Annular Modes.
!     J. Atmos. Sci., 61, 264-280.  (equation A.6)
!     The forcing at each wavenumber follows a Markov process with decorrelation
!     timescale 'decay_time'. The spectral coefficients are projected to physical
!     space and multiplied by a Gaussian latitude (and optionally longitude) mask
!     before being added to the vorticity tendency.
!
! (2) White-in-time stochastic forcing  [do_white_forcing = .true.]
!     After:
!     Scott, R. K. and Polvani, L. M., 2007: Forced-dissipative shallow-water
!     turbulence on the sphere and the atmospheric circulation of the giant
!     planets. J. Atmos. Sci., 64, 3158-3176.
!
!     Each timestep, fresh independent random complex spectral coefficients
!     with FIXED AMPLITUDE and UNIFORMLY RANDOM PHASE are drawn for all
!     modes in the forcing band (delta-correlated in time). This exactly
!     matches S&P's forcing structure: f_mn = amplitude * exp(i*theta),
!     theta ~ U[0, 2*pi).
!
!     The forcing amplitude is set via epsilon_0 (m^2/s^3), the kinetic
!     energy injection rate per unit mass. 'amplitude' is computed
!     internally from epsilon_0 and must not be set in the namelist.
!
!     The formula follows S&P equation (1) exactly, converted to dimensional
!     vorticity tendency units and corrected for the GFDL leapfrog scheme
!     (atmosphere.F90 passes dt_real to stirring_init but the leapfrog step
!     is 2*dt_real, introducing a factor of sqrt(2)):
!
!       amplitude = (1/a) * sqrt( 4 * l_f*(l_f+1) * epsilon_0 /
!                                     ( (2*l_f+1) * dt * delta_n ) )
!
!     where:
!       l_f     = (n_total_forcing_min + n_total_forcing_max) / 2
!                 (central total wavenumber of the forcing band)
!       delta_n = n_total_forcing_max - n_total_forcing_min - 1
!                 (number of forced total wavenumbers)
!       dt      = dt_real as passed to stirring_init (= dt_atmos from namelist)
!       a       = radius  [from constants_mod]
!
!     The factor of 4 = 2 (S&P) * 2 (leapfrog correction).
!     S&P's original factor of 2 accounts for the two complex degrees of
!     freedom (real and imaginary parts) per mode. The additional factor of 2
!     corrects for the GFDL leapfrog, verified empirically.
!
!     epsilon_0 has units of m^2/s^3 (dimensional). To convert from the
!     non-dimensional S&P value eps* (non-dimensionalised by radius^2*omega^3):
!       epsilon_0 = eps* * radius^2 * omega^3
!     where omega is the planetary rotation rate (rad/s).
!
!     For isotropic S&P-style forcing also set do_localize=.false. and
!     zonal_forcing_min=0.

! Stirring is not part of barotropic_physics because barotropic_physics appears to be intended for
! operations that are done completely in grid space. Stirring is computed partly in spectral space.

use    constants_mod, only: pi, radius

use time_manager_mod, only: time_type

use          fms_mod, only: open_namelist_file, check_nml_error, close_file, write_version_number, &
                            stdlog, mpp_pe, mpp_root_pe, file_exist, read_data, write_data, error_mesg, FATAL

use   transforms_mod, only: get_spec_domain, get_grid_domain, trans_spherical_to_grid, trans_grid_to_spherical, &
                            grid_domain, get_lon_max, get_lat_max, get_deg_lon, get_deg_lat, get_grid_boundaries, &
                            get_num_fourier, get_num_spherical, spectral_domain

use diag_manager_mod, only: diag_axis_init, register_static_field, register_diag_field, send_data

implicit none
private

integer :: ms,me,ns,ne,is,ie,js,je
integer :: id_str_amp, id_g_stir_sqr, id_stir
logical :: used
logical, allocatable, dimension(:,:) :: wave_mask   ! wave_mask(m,n) = .true. if spherical wave (m,n) is to be excited
complex, allocatable, dimension(:,:) :: s_stir      ! stirring. Saved from one time step to the next
real,    allocatable, dimension(:,:) :: localize    ! localizes the stirring
real,    allocatable, dimension(:,:) :: g_stir_sqr  ! time mean of g_stir**2 over entire integration
integer, allocatable, dimension(:)   :: seed        ! random number seed
real ::  astir, bstir
integer :: num_steps, num_fourier, num_spherical, nseed

logical :: module_is_initialized = .false.

character(len=128) :: version = '$Id: stirring.F90,v 17.0 2009/07/21 03:00:25 fms Exp $'
character(len=128) :: tagname = '$Name: siena_201207 $'

public :: stirring_init, stirring, stirring_end

! -----------------------------------------------------------------------
! Namelist parameters
! -----------------------------------------------------------------------

! OU forcing (do_white_forcing=.false.) parameters:
real :: decay_time = 2*86400  ! OU decorrelation timescale (s)
real :: amplitude  = 0.0      ! vorticity tendency amplitude (s^-2).
                               ! For white-in-time forcing, leave at 0.0 and
                               ! set epsilon_0 instead; amplitude is then computed internally.
real :: lat0   = 45.          ! centre latitude of Gaussian mask (degrees)
real :: widthy = 12.          ! meridional half-width of Gaussian mask (degrees)
real :: lon0   = 180.         ! centre longitude of zonal structure (degrees)
real :: widthx = 45.          ! zonal half-width of Gaussian structure (degrees)
real :: B      = 0.0          ! relative amplitude of zonal structure (0 = zonally symmetric)
real :: C      = 1.0          ! (unused, retained for compatibility)

logical :: do_localize = .true.  ! .true.: apply Gaussian lat/lon mask in physical space (Vallis)
                                  ! .false.: isotropic forcing everywhere (S&P)

! Forcing wavenumber band. Forced total wavenumbers L = m+n satisfy:
!   n_total_forcing_min < L < n_total_forcing_max   AND   m > zonal_forcing_min
integer :: n_total_forcing_max = 15  ! total wavenumbers strictly less than this are forced
integer :: n_total_forcing_min = 9   ! total wavenumbers strictly greater than this are forced
integer :: zonal_forcing_min   = 3   ! zonal wavenumbers strictly greater than this are forced
                                      ! Set to 0 for isotropic S&P forcing (all m >= 1).

! White-in-time forcing option (Scott & Polvani 2007):
logical :: do_white_forcing = .false.

! Energy injection rate for white-in-time forcing (m^2/s^3, per unit mass).
! When do_white_forcing=.true., set epsilon_0 > 0 and leave amplitude = 0.0.
! amplitude is then computed internally from epsilon_0 (see module header).
real :: epsilon_0 = 0.0

! If >= 0, deterministically seeds the stirring random number generator from this value, giving
! reproducible stirring on a fresh (non-restart) run. Ignored when restarting, since the restart
! file's own saved seed is used instead. Default -1 leaves the runtime's own default (non-deterministic
! on a fresh run) seed in place, matching the pre-existing behaviour.
integer :: fixed_random_seed = -1

namelist / stirring_nml / decay_time, amplitude, lat0, lon0, widthy, widthx, B, do_localize, &
                          n_total_forcing_max, n_total_forcing_min, zonal_forcing_min, &
                          do_white_forcing, epsilon_0, fixed_random_seed

contains

!================================================================================================================================
subroutine stirring_init(dt, Time, id_lon, id_lat, id_lonb, id_latb)
real,             intent(in) :: dt         ! = dt_atmos (s), as passed from atmosphere.F90
type(time_type),  intent(in) :: Time
integer,          intent(in) :: id_lon, id_lat, id_lonb, id_latb

real    :: xx, kk
real    :: l_f, delta_n  ! S&P central wavenumber and band width; used in amplitude formula
integer :: i, j, m, n, L, ierr, io, unit, lon_max, lat_max
real,    allocatable, dimension(:)   :: ampx, ampy, lon, lat, lonb, latb
real,    allocatable, dimension(:,:) :: real_part, imag_part

if(module_is_initialized) return

call write_version_number (version, tagname)

if (file_exist('input.nml')) then
  unit = open_namelist_file ()
  ierr = 1
  do while (ierr /= 0)
    read  (unit, nml=stirring_nml, iostat=io, end=10)
    ierr = check_nml_error (io, 'stirring_nml')
  enddo
  10 call close_file (unit)
endif
if(mpp_pe() == mpp_root_pe()) write(stdlog(), nml=stirring_nml)

! -----------------------------------------------------------------------
! Compute amplitude from epsilon_0 (white-in-time forcing only).
!
! This block must come BEFORE the 'if(amplitude==0.0) return' guard so
! that epsilon_0 alone is sufficient to activate the forcing.
!
! Exact S&P formula (their eq. 1), converted to dimensional vorticity
! tendency and corrected for the GFDL leapfrog (see module header):
!
!   amplitude = (1/a) * sqrt( 4 * l_f*(l_f+1) * epsilon_0 /
!                                 ( (2*l_f+1) * dt * delta_n ) )
!
! where l_f = (n_total_forcing_min + n_total_forcing_max) / 2
!       delta_n = n_total_forcing_max - n_total_forcing_min - 1
! -----------------------------------------------------------------------
if (do_white_forcing .and. epsilon_0 > 0.0) then
  if (amplitude /= 0.0) then
    call error_mesg('stirring_init', &
      'Set either epsilon_0 OR amplitude, not both. Set amplitude=0.0 when using epsilon_0.', FATAL)
  endif
  delta_n = real(n_total_forcing_max - n_total_forcing_min - 1)
  if (delta_n <= 0.0) then
    call error_mesg('stirring_init', &
      'Forcing band is empty: n_total_forcing_max must be > n_total_forcing_min + 1.', FATAL)
  endif
  l_f = real(n_total_forcing_min + n_total_forcing_max) / 2.0
  amplitude = (1.0/radius) * sqrt(4.0 * l_f*(l_f+1.0) * epsilon_0 &
                                   / ((2.0*l_f+1.0) * dt * delta_n))
  if(mpp_pe() == mpp_root_pe()) then
    write(stdlog(),'(a,es13.5,a)') 'stirring_mod: epsilon_0 = ', epsilon_0, ' m^2/s^3'
    write(stdlog(),'(a,f8.1)')     'stirring_mod: l_f       = ', l_f
    write(stdlog(),'(a,f8.1)')     'stirring_mod: delta_n   = ', delta_n
    write(stdlog(),'(a,es13.5,a)') 'stirring_mod: amplitude = ', amplitude, ' s^-2  (computed from epsilon_0, S&P eq.1)'
  endif
endif

call get_lon_max(lon_max)
call get_lat_max(lat_max)

allocate(lon (lon_max  )) ; lon  = 0.0
allocate(lat (lat_max  )) ; lat  = 0.0
allocate(lonb(lon_max+1)) ; lonb = 0.0
allocate(latb(lat_max+1)) ; latb = 0.0

call get_deg_lon(lon)
call get_deg_lat(lat)

module_is_initialized = .true.
if(amplitude == 0.0) return  ! stirring does nothing more unless amplitude is non-zero

call get_spec_domain(ms,me,ns,ne)
call get_grid_domain(is,ie,js,je)
call get_num_fourier(num_fourier)
call get_num_spherical(num_spherical)

allocate(wave_mask(ms:me,ns:ne)) ; wave_mask   = .false.
allocate(s_stir  (ms:me,ns:ne)) ; s_stir      = cmplx(0.0,0.0)
allocate(ampx    (is:ie))        ; ampx        = 0.0
allocate(ampy    (js:je))        ; ampy        = 0.0
allocate(localize(is:ie,js:je))  ; localize    = 0.0
allocate(g_stir_sqr(is:ie,js:je)); g_stir_sqr  = 0.0

! wave_mask(m,n) = .true. when:
!   (m+n) > n_total_forcing_min  .AND.  (m+n) < n_total_forcing_max  .AND.  m > zonal_forcing_min
! For S&P isotropic forcing: set zonal_forcing_min=0 so all m >= 1 are forced.
do m = (zonal_forcing_min+1), (n_total_forcing_max-1)
  if(m >= ms .and. m <= me) then
    do n = (n_total_forcing_min+1)-m, (n_total_forcing_max-1)-m
      if(n >= ns .and. n <= ne) then
        wave_mask(m,n) = .true.
      endif
    enddo
  endif
enddo

! OU coefficients (unused for white-in-time forcing; set to safe values).
if (do_white_forcing) then
  astir = 1.0
  bstir = 0.0
  if(mpp_pe() == mpp_root_pe()) then
    write(stdlog(),'(a)') &
      'stirring_mod: do_white_forcing=T -> white-in-time (S&P) forcing active.'
    write(stdlog(),'(a)') &
      'stirring_mod: fixed amplitude, random phase: |f_mn|=amplitude, phase~U[0,2*pi].'
    write(stdlog(),'(a)') &
      'stirring_mod: decay_time is ignored.'
  endif
else
  astir = sqrt(1.0 - exp(-2.0*dt/decay_time))
  bstir = exp(-dt/decay_time)
endif

do i = is, ie
  xx = lon(i) - lon0
  xx = xx - 360.*nint(xx/360.)  ! wrap to [-180, +180]
  ampx(i) = 1.0 + B*exp(-0.5*(xx/widthx)**2)
enddo
do j = js, je
  ampy(j) = exp(-0.5*((lat(j)-lat0)/widthy)**2)
enddo
if (do_localize) then
  do j = js, je
    do i = is, ie
      localize(i,j) = ampx(i)*ampy(j)
    enddo
  enddo
else
  localize = 1.0
endif

deallocate(ampx, ampy)

num_steps = 0
id_g_stir_sqr = register_static_field('stirring_mod', 'stirring_sqr', (/id_lon,id_lat/), &
                                       'stirring squared', '1/sec^4')
id_str_amp    = register_static_field('stirring_mod', 'stirring_amp', (/id_lon,id_lat/), &
                                       'amplitude of stirring', 'none')
id_stir       = register_diag_field  ('stirring_mod', 'stirring',     (/id_lon,id_lat/), &
                                       Time, 'stirring', '1/sec^2')
used = send_data(id_str_amp, amplitude*localize)

call random_seed(size=nseed)
allocate(seed(nseed))

if(file_exist('INPUT/stirring.res.nc')) then
  allocate(real_part(ms:me,ns:ne), imag_part(ms:me,ns:ne))
  call read_data('INPUT/stirring.res.nc', 'stir_real', real_part, spectral_domain)
  call read_data('INPUT/stirring.res.nc', 'stir_imag', imag_part, spectral_domain)
  do n = ns, ne
    do m = ms, me
      s_stir(m,n) = cmplx(real_part(m,n), imag_part(m,n))
    end do
  end do
  deallocate(real_part, imag_part)
  call read_data('INPUT/stirring.res.nc', 'ran_nmbr_seed', seed, no_domain=.true.)
  call random_seed(put=seed)
else if(fixed_random_seed >= 0) then
  seed = fixed_random_seed + (/ (i, i=1,nseed) /)
  call random_seed(put=seed)
endif
! Note: for white-in-time forcing, the s_stir read from restart is overwritten on the
! first call to stirring(), so it is irrelevant. The random seed IS restored, ensuring
! reproducible restart runs.

end subroutine stirring_init
!================================================================================================================================
subroutine stirring(Time, dt_vors)
type(time_type),                     intent(in)    :: Time
complex, dimension(ms:me,ns:ne),     intent(inout) :: dt_vors

real,    dimension(is:ie,js:je)             :: g_stir
complex, dimension(ms:me,ns:ne)             :: new_stirring
real,    dimension(0:num_fourier,0:num_spherical,2) :: ran_nmbrs
real    :: theta
integer :: m, n

if(.not.module_is_initialized) then
  call error_mesg('stirring', 'stirring_init has not been called', FATAL)
end if

if(amplitude == 0.0) return  ! stirring does nothing unless amplitude is non-zero

call random_number(ran_nmbrs)

! -----------------------------------------------------------------------
! Generate spectral forcing coefficients for all forced modes.
!
! Vallis OU (do_white_forcing=.false.):
!   new_stirring = amplitude * astir * cmplx(U[-1,1], U[-1,1])
!   astir = sqrt(1 - exp(-2*dt/decay_time)) gives the correct OU steady-state
!   variance after blending with the previous s_stir below.
!
! White-in-time (do_white_forcing=.true.) [Scott & Polvani 2007]:
!   new_stirring = amplitude * exp(i*theta), theta = 2*pi * U[0,1)
!   Fixed magnitude, uniformly random phase: exactly S&P's formulation.
!   Only ran_nmbrs(:,:,1) is used (phase); ran_nmbrs(:,:,2) is unused.
! -----------------------------------------------------------------------
do n = ns, ne
  do m = ms, me
    if(wave_mask(m,n)) then
      if(do_white_forcing) then
        theta = 2.0 * pi * ran_nmbrs(m,n,1)
        new_stirring(m,n) = amplitude * cmplx(cos(theta), sin(theta))
      else
        new_stirring(m,n) = amplitude * astir * cmplx(2*ran_nmbrs(m,n,1)-1, 2*ran_nmbrs(m,n,2)-1)
      endif
    else
      new_stirring(m,n) = cmplx(0.0, 0.0)
    endif
  enddo
enddo

! -----------------------------------------------------------------------
! Apply physical-space spatial localization (Vallis-style only).
! For S&P isotropic forcing: set do_localize=.false. to skip the transform
! round-trip entirely (saves 2 spectral transforms per timestep).
! When do_localize=.true., the physical-space grid-space Gaussian mask is
! applied, and the global-mean artefact introduced by this is zeroed below.
! -----------------------------------------------------------------------
if (do_localize) then
  call trans_spherical_to_grid(new_stirring, g_stir)
  g_stir = localize * g_stir
  call trans_grid_to_spherical(g_stir, new_stirring)
  if(ms == 0 .and. ns == 0) then
    new_stirring(0,0) = cmplx(0.0,0.0) ! zero the global mean artefact
  endif
endif

! -----------------------------------------------------------------------
! Update s_stir and add to the vorticity tendency.
!
! Vallis OU: s_stir = bstir*s_stir + new_stirring  (Vallis et al. 2004, eq. A.6)
!   bstir = exp(-dt/decay_time) gives exponentially autocorrelated forcing.
!
! White-in-time: s_stir = new_stirring
!   No memory: forcing is delta-correlated in time.
! -----------------------------------------------------------------------
if (do_white_forcing) then
  s_stir = new_stirring
else
  s_stir = bstir*s_stir + new_stirring  ! equation A.6, Vallis et al. 2004
                                         ! DOI:10.1175/1520-0469(2004)061<0264:AMASDM>2.0.CO;2
endif

dt_vors = dt_vors + s_stir
call trans_spherical_to_grid(s_stir, g_stir)
g_stir_sqr = g_stir_sqr + g_stir*g_stir
num_steps  = num_steps + 1
used = send_data(id_stir, g_stir, Time)

end subroutine stirring
!================================================================================================================================
subroutine stirring_end

if(.not.module_is_initialized) return
if(amplitude == 0.0) return  ! stirring does nothing unless amplitude is non-zero

g_stir_sqr = g_stir_sqr / num_steps
used = send_data(id_g_stir_sqr, g_stir_sqr)

call write_data('RESTART/stirring.res.nc', 'stir_real',  real(s_stir), spectral_domain)
call write_data('RESTART/stirring.res.nc', 'stir_imag', aimag(s_stir), spectral_domain)
call random_seed(get=seed)
call write_data('RESTART/stirring.res.nc', 'ran_nmbr_seed', seed, no_domain=.true.)

deallocate(wave_mask, s_stir, localize, g_stir_sqr)
module_is_initialized = .false.

end subroutine stirring_end
!================================================================================================================================

end module stirring_mod