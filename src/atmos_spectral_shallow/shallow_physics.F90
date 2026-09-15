module shallow_physics_mod

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

use               fms_mod, only: open_namelist_file,   &
                                 open_restart_file,    &
                                 file_exist,           &
                                 check_nml_error,      &
                                 error_mesg,           &
                                 FATAL, WARNING,       &
                                 write_version_number, &
                                 mpp_pe,               &
                                 mpp_root_pe,          &
                                 fms_init, fms_end,    &
                                 read_data,            &
                                 write_data,           &
                                 set_domain,           &
                                 close_file,           &
                                 stdlog

use         transforms_mod, only: get_sin_lat, get_cos_lat,  &
                                  get_deg_lon, get_deg_lat,  &
                                  get_wts_lat, &
                                  get_grid_domain, get_spec_domain, &
                                  grid_domain, area_weighted_global_mean

use       time_manager_mod, only: time_type

use         constants_mod, only: omega

use   shallow_dynamics_mod, only: get_u_deep_mag

!========================================================================
implicit none
private
!========================================================================

public :: shallow_physics_init,    &
          shallow_physics,         &
          shallow_physics_end,     &
          phys_type


! version information 
!========================================================================
character(len=128) :: version = '$Id: shallow_physics.F90,v 10.0 2003/10/24 22:01:02 fms Exp $'
character(len=128) :: tagname = '$Name: siena_201207 $'
!========================================================================

type phys_type
   real, pointer, dimension(:,:)   :: empty=>NULL()
   real, pointer, dimension(:,:)   :: h_eq=>NULL()
   real, pointer, dimension(:,:)   :: du_dt_mass=>NULL()
   real, pointer, dimension(:,:)   :: dv_dt_mass=>NULL()
   real, pointer, dimension(:,:)   :: du_dt_drag=>NULL()
   real, pointer, dimension(:,:)   :: dv_dt_drag=>NULL()
   real, pointer, dimension(:,:)   :: u_deep=>NULL()
end type

logical :: module_is_initialized = .false.

integer :: is, ie, js, je

integer :: pe
logical :: root

real, allocatable, dimension(:) :: rad_lat, deg_lat, deg_lon, &
         sin_lat, cos_lat, wts_lat, coriolis

real, allocatable, target, dimension(:,:) :: h_eq, du_dt_mass, dv_dt_mass, &
         du_dt_drag, dv_dt_drag, u_deep

real    :: kappa_m, kappa_t

real    :: h_warn                 ! soft floor for the low h warning, 0 disables
logical :: warned_low_h = .false.



! namelist 
!========================================================================

real    :: fric_damp_time  = -20.0
real    :: therm_damp_time = -10.0
real    :: del_h           = 0.0
real    :: h_0             = 3.e04
real    :: h_amp           = 2.e04
real    :: h_lon           =  90.0
real    :: h_lat           =  25.0
real    :: h_width         =  15.0
real    :: h_itcz          = 1.e05
real    :: itcz_width      =  4.0

! h_eq_option: 'legacy' (gaussian blob + itcz), 'showman_polvani'
! (h_0 + del_h*h_0*cos(lat)*cos(lon-h_lon)), or 'perez_becker' (as
! showman_polvani on the dayside, flat on the nightside). For the two tidally
! locked options h_lon is the substellar longitude (set it to 0.0; the default
! of 90.0 is the legacy blob centre), the substellar latitude is 0, and del_h
! is the day-night contrast as a fraction of h_0, not a geopotential.
character(len=64) :: h_eq_option = 'legacy'

logical :: do_zero_mean_h_eq = .true.   ! shift h_eq so its global mean is h_0
logical :: do_mass_exchange  = .false.  ! -Q*u/h and -Q*v/h tendencies where Q > 0

! deep jets. u_deep ~ u_deep_amp*E(lat)*cos(u_deep_n*lat), with lat in radians and
! E a gaussian taper, disabled if its width is <= 0
logical :: do_deep_jet_force = .false.  ! f*u_deep body force on the vcomp tendency
logical :: do_deep_jet_mass  = .false.  ! mass exchange injects u_deep, not zero
real    :: u_deep_amp   = 200.0
real    :: u_deep_n     =   8.0
real    :: u_deep_width =  20.0         ! degrees

namelist /shallow_physics_nml/ fric_damp_time, therm_damp_time, del_h, h_0, &
                               h_amp, h_lon, h_lat, h_width, &
                               itcz_width, h_itcz, h_eq_option, &
                               do_zero_mean_h_eq, do_mass_exchange, &
                               do_deep_jet_force, do_deep_jet_mass, &
                               u_deep_amp, u_deep_n, u_deep_width
!========================================================================

contains

!========================================================================

subroutine shallow_physics_init(Phys) 

type(phys_type), intent(inout) :: Phys

integer :: i, j, unit, ierr, io

real :: xx, yy, dd, coszen, h_eq_mean, c_deep

logical :: dayside_only

real, allocatable, dimension(:)   :: taper
real, allocatable, dimension(:,:) :: am_num, am_den

call write_version_number(version, tagname)

pe = mpp_pe()
root = (pe == mpp_root_pe())

! read the namelist

if (file_exist('input.nml')) then
  unit = open_namelist_file ()
  ierr=1
  do while (ierr /= 0)
    read  (unit, nml=shallow_physics_nml, iostat=io, end=10)
    ierr = check_nml_error (io, 'shallow_physics_nml')
  enddo
  10 call close_file (unit)
endif

if(fric_damp_time  < 0.0)  fric_damp_time = -  fric_damp_time*86400
if(therm_damp_time < 0.0) therm_damp_time = - therm_damp_time*86400

kappa_m = 0.0
kappa_t = 0.0
if( fric_damp_time .ne. 0.0) kappa_m = 1./fric_damp_time
if(therm_damp_time .ne. 0.0) kappa_t = 1./therm_damp_time

call get_grid_domain(is,ie,js,je)

allocate ( rad_lat      (js:je) )
allocate ( deg_lat      (js:je) )
allocate ( sin_lat      (js:je) )
allocate ( cos_lat      (js:je) )
allocate ( wts_lat      (js:je) )
allocate ( coriolis     (js:je) )
allocate ( deg_lon      (is:ie) )
allocate ( h_eq   (is:ie,js:je) )
allocate ( u_deep (is:ie,js:je) )
allocate ( du_dt_mass (is:ie,js:je) ) ; du_dt_mass = 0.0
allocate ( dv_dt_mass (is:ie,js:je) ) ; dv_dt_mass = 0.0
allocate ( du_dt_drag (is:ie,js:je) ) ; du_dt_drag = 0.0
allocate ( dv_dt_drag (is:ie,js:je) ) ; dv_dt_drag = 0.0

call get_wts_lat(wts_lat)
call get_deg_lat(deg_lat)
call get_deg_lon(deg_lon)
rad_lat = deg_lat*atan(1.)/45.
sin_lat = sin(rad_lat)
cos_lat = cos(rad_lat)
coriolis = 2*omega*sin_lat


select case (trim(h_eq_option))

case ('legacy')

  do j = js, je
    do i = is, ie
       xx = (deg_lon(i) - h_lon)/(h_width*2.0)
       yy = (deg_lat(j) - h_lat)/h_width
       dd =  xx*xx + yy*yy
       h_eq(i,j) = h_0 + h_amp*max(1.e-10, exp(-dd))
    end do
  end do

  do j = js, je
    yy = deg_lat(j)/itcz_width
    dd = yy*yy
    h_eq(:,j) = h_eq(:,j) + h_itcz*exp(-dd)
  end do

  h_warn = 0.0

case ('showman_polvani', 'perez_becker')

  dayside_only = (trim(h_eq_option) == 'perez_becker')

  do j = js, je
    do i = is, ie
      ! cosine of the angle from the substellar point at (lat=0, lon=h_lon)
      coszen = cos_lat(j)*cos((deg_lon(i) - h_lon)*atan(1.)/45.)
      if(dayside_only) coszen = max(coszen, 0.0)
      h_eq(i,j) = h_0 + del_h*h_0*coszen
    end do
  end do

  ! a uniform shift, so the day-night contrast is del_h*h_0 either way
  if(do_zero_mean_h_eq) then
    h_eq_mean = area_weighted_global_mean(h_eq)
    h_eq = h_eq - (h_eq_mean - h_0)
  endif

  h_warn = 0.05*h_0

case default

  call error_mesg('shallow_physics_init', &
                  'unrecognised h_eq_option: '//trim(h_eq_option), FATAL)

end select

if(minval(h_eq) <= 0.0) then
  call error_mesg('shallow_physics_init', 'h_eq is not positive everywhere', FATAL)
endif

allocate ( taper       (js:je) )
allocate ( am_num (is:ie,js:je) )
allocate ( am_den (is:ie,js:je) )

do j = js, je
  taper(j) = 1.0
  if(u_deep_width > 0.0) taper(j) = exp(-deg_lat(j)**2/(2.*u_deep_width**2))
  u_deep(:,j) = taper(j)*cos(u_deep_n*rad_lat(j))
  am_num(:,j) = u_deep(:,j)*cos_lat(j)
  am_den(:,j) = taper(j)*cos_lat(j)
end do

! area_weighted_global_mean carries one cos(lat) in its weights, so the explicit
! cos_lat above gives the cos(lat)^2 angular momentum weighting, and the sum runs
! over all processors
c_deep = area_weighted_global_mean(am_num)/area_weighted_global_mean(am_den)

! subtracting c_deep leaves the deep jets with zero integrated angular momentum
do j = js, je
  u_deep(:,j) = u_deep_amp*(u_deep(:,j) - taper(j)*c_deep)
end do

! u_deep is unused, so set to zero for diagnostics
if(.not.do_deep_jet_force .and. .not.do_deep_jet_mass) u_deep = 0.0

deallocate ( taper, am_num, am_den )

if(do_deep_jet_force .and. get_u_deep_mag() /= 0.0) then
  call error_mesg('shallow_physics_init', &
                  'do_deep_jet_force and shallow_dynamics_nml u_deep_mag are both '// &
                  'active, which applies the deep jet force twice', FATAL)
endif

if(do_deep_jet_mass .and. .not.do_mass_exchange) then
  call error_mesg('shallow_physics_init', &
                  'do_deep_jet_mass does nothing unless do_mass_exchange is true', FATAL)
endif

Phys%h_eq       => h_eq
Phys%du_dt_mass => du_dt_mass
Phys%dv_dt_mass => dv_dt_mass
Phys%du_dt_drag => du_dt_drag
Phys%dv_dt_drag => dv_dt_drag
Phys%u_deep     => u_deep

!if(file_exist('INPUT/shallow_physics.res')) then
!  unit = open_restart_file(file='INPUT/shallow_physics.res',action='read')
!  call set_domain(grid_domain)
!  call close_file(unit)
!else

!endif

module_is_initialized = .true.

return
end subroutine shallow_physics_init

!=======================================================================

subroutine shallow_physics(Time, dt_ug, dt_vg, dt_hg, ug, vg, hg,   &
                             delta_t, previous, current, Phys)

real, intent(inout),  dimension(is:ie, js:je)    :: dt_ug, dt_vg, dt_hg
real, intent(in)   ,  dimension(is:ie, js:je, 2) :: ug, vg, hg

real   , intent(in)  :: delta_t
integer, intent(in)  :: previous, current

type(time_type), intent(in)    :: Time
type(phys_type), intent(inout) :: Phys

real, dimension(is:ie, js:je) :: q_mass, q_rate

real    :: h_min_local
integer :: j

h_min_local = minval(hg(:,:,previous))

if(h_min_local <= 0.0) then
  call error_mesg('shallow_physics', 'layer geopotential is not positive', FATAL)
endif

if(.not.warned_low_h .and. h_min_local < h_warn) then
  call error_mesg('shallow_physics', &
                  'layer geopotential has fallen below 0.05*h_0', WARNING)
  warned_low_h = .true.
endif

! mass source, reused by the mass exchange term below
q_mass = kappa_t*(h_eq - hg(:,:,previous))

du_dt_drag = -kappa_m*ug(:,:,previous)
dv_dt_drag = -kappa_m*vg(:,:,previous)

dt_ug = dt_ug + du_dt_drag
dt_vg = dt_vg + dv_dt_drag
dt_hg = dt_hg + q_mass

! showman and polvani mass exchange: injected mass carries no momentum, or the
! deep jet momentum if do_deep_jet_mass
if(do_mass_exchange) then
  q_rate     = max(q_mass, 0.0)/hg(:,:,previous)
  du_dt_mass = -q_rate*ug(:,:,previous)
  if(do_deep_jet_mass) du_dt_mass = du_dt_mass + q_rate*u_deep
  dv_dt_mass = -q_rate*vg(:,:,previous)
  dt_ug      = dt_ug + du_dt_mass
  dt_vg      = dt_vg + dv_dt_mass
endif

! deep jet pressure gradient felt by the layer, f*k x u_deep
if(do_deep_jet_force) then
  do j = js, je
    dt_vg(:,j) = dt_vg(:,j) + coriolis(j)*u_deep(:,j)
  end do
endif


return
end subroutine shallow_physics

!======================================================================

subroutine shallow_physics_end(Phys)

type(phys_type), intent(in) :: Phys

integer :: unit

if(.not.module_is_initialized) then
  call error_mesg('shallow_physics_end','physics has not been initialized ', FATAL)
endif

!unit = open_restart_file(file='RESTART/shallow_physics.res', action='write')

!call set_domain(grid_domain)

!call close_file(unit)

module_is_initialized = .false.

return
end subroutine shallow_physics_end

!======================================================================

end module shallow_physics_mod
