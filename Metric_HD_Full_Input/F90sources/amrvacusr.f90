!=============================================================================
! amrvacusr.t
!=============================================================================
!  INCLUDE:amrvacnul/speciallog.t
!  INCLUDE:amrvacnul/specialbound.t
!  INCLUDE:amrvacnul/specialsource.t
!================================================================================
!
!    BHAC (The Black Hole Accretion Code) solves the equations of
!    general relativistic magnetohydrodynamics and other hyperbolic systems
!    in curved spacetimes.
!
!    Copyright (C) 2019 Oliver Porth, Hector Olivares, Yosuke Mizuno, Ziri Younsi,
!    Luciano Rezzolla, Elias Most, Bart Ripperda and Fabio Bacchini
!
!    This file is part of BHAC.
!
!    BHAC is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    BHAC is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with BHAC.  If not, see <https://www.gnu.org/licenses/>.
!
!================================================================================

!=============================================================================
subroutine specialsource_impl(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iwmin,iwmax,qtC,&
   wCT,qt,w,x)

use mod_amrvacdef

integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
    ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, iwmin,iwmax
double precision, intent(in) :: qdt, qtC, qt, x(ixImin1:ixImax1,&
   ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw), wCT(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nw)
!-----------------------------------------------------------------------------

end subroutine specialsource_impl
!=============================================================================
subroutine getdt_impl(w,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
   ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,dtnew,dx1,dx2,dx3,x)

use mod_amrvacdef

integer, intent(in) :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,ixmin1,&
   ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
double precision, intent(in) :: dx1,dx2,dx3, x(ixGmin1:ixGmax1,&
   ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:ndim)
! note that depending on strictsmall etc, w values may change
double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
   ixGmin3:ixGmax3,1:nw), dtnew
!-----------------------------------------------------------------------------
dtnew=bigdouble

end subroutine getdt_impl
!=============================================================================
!  INCLUDE:amrvacnul/usrflags.t
!  INCLUDE:amrvacnul/correctaux_usr.t
  !=============================================================================
  subroutine initglobaldata_usr

    use mod_amrvacdef
    use mod_oneblock
    !-----------------------------------------------------------------------------

    eqpar(gamma_)  = 5.0d0/3.0d0 ! Adiabatic index
    
    call read_oneblock("it_18432.blk")

!    call mpistop('stopping after reading')

  end subroutine initglobaldata_usr
  !=============================================================================
  subroutine initonegrid_usr(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,s)

    ! initialize one grid within ixO^L

    use mod_amrvacdef
    use mod_oneblock

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
    type(state)         :: s
    ! .. local ..
    integer             :: ix1,ix2,ix3
    !-----------------------------------------------------------------------------
    associate(x=>s%x%x,w=>s%w%w)

       w = 0.0d0
      
       
       
do ix1=ixOmin1, ixOmax1 
do ix2=ixOmin2, ixOmax2 
do ix3=ixOmin3, ixOmax3 
    call interpolate_oneblock( x(ix1,ix2,ix3,:) , rho_, w(ix1,ix2,ix3, rho_) )
    call interpolate_oneblock( x(ix1,ix2,ix3,:) , pp_, w(ix1,ix2,ix3, pp_) )

    call interpolate_oneblock( x(ix1,ix2,ix3,:) , u1_, w(ix1,ix2,ix3, u1_) )
    call interpolate_oneblock( x(ix1,ix2,ix3,:) , u2_, w(ix1,ix2,ix3, u2_) )
    call interpolate_oneblock( x(ix1,ix2,ix3,:) , u3_, w(ix1,ix2,ix3, u3_) )

enddo
enddo
enddo

       
      
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,s_) &
         = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         pp_) * w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)&
         **(-eqpar(gamma_))
      

    call conserve(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,patchfalse)
 
  end associate
end subroutine initonegrid_usr
!=============================================================================
subroutine initvecpot_usr(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
    ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, xC, A, idir)

  ! initialize the vectorpotential on the corners
  ! used by b_from_vectorpotential()


  use mod_amrvacdef

  integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3, ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, idir
  double precision, intent(in)       :: xC(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:ndim)
  double precision, intent(out)      :: A(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3)
  ! .. local ..
  !-----------------------------------------------------------------------------

end subroutine initvecpot_usr
!=============================================================================
subroutine specialvar_output(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,nwmax,w,s,normconv)

! this subroutine can be used in convert, to add auxiliary variables to the
! converted output file, for further analysis using tecplot, paraview, ....
! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
!
! the array normconv can be filled in the (nw+1:nw+nwauxio) range with 
! corresponding normalization values (default value 1)

use mod_amrvacdef

integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,nwmax
double precision                   :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nwmax)
type(state)                        :: s
double precision                   :: normconv(0:nwmax)
!-----------------------------------------------------------------------------
associate(x=>s%x%x)


! Reduce output array size, +1 was added for eventual pointdata output
call get_divb(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1+1,&
   ixOmin2+1,ixOmin3+1,ixOmax1-1,ixOmax2-1,ixOmax3-1,w(ixImin1:ixImax1,&
   ixImin2:ixImax2,ixImin3:ixImax3,1:nw),w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,nw+1))


if (nwmax >= nw+2) then
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+2) &
      = myM%alpha(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
end if

if (nwmax >= nw+3) then
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+3) &
      = myM%g(1,1)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
end if

if (nwmax >= nw+4) then
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+4) &
      = myM%g(1,2)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
end if

if (nwmax >= nw+5) then
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+5) &
      = myM%g(1,3)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
end if

if (nwmax >= nw+6) then
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+6) &
      = myM%g(2,2)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
end if

if (nwmax >= nw+7) then
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+7) &
      = myM%g(2,3)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
end if

if (nwmax >= nw+8) then
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+8) &
      = myM%g(3,3)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
end if

if (nwmax >= nw+9) then
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+9) &
      = myM%beta(1)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
end if

if (nwmax >= nw+10) then
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+10) &
      = myM%beta(2)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
end if

if (nwmax >= nw+11) then
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+11) &
      = myM%beta(3)%elem(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
end if

end associate
end subroutine specialvar_output
!=============================================================================
subroutine specialvarnames_output

! newly added variables to be concatenated with the primnames/wnames string

use mod_amrvacdef
!-----------------------------------------------------------------------------

primnames= TRIM(primnames)//' '//'divB'
wnames=TRIM(wnames)//' '//'divB'

primnames= TRIM(primnames)//' '//'alpha'
wnames=TRIM(wnames)//' '//'alpha'

primnames= TRIM(primnames)//' '//'g11'
wnames=TRIM(wnames)//' '//'g11'

primnames= TRIM(primnames)//' '//'g12'
wnames=TRIM(wnames)//' '//'g12'

primnames= TRIM(primnames)//' '//'g13'
wnames=TRIM(wnames)//' '//'g13'

primnames= TRIM(primnames)//' '//'g22'
wnames=TRIM(wnames)//' '//'g22'

primnames= TRIM(primnames)//' '//'g23'
wnames=TRIM(wnames)//' '//'g23'

primnames= TRIM(primnames)//' '//'g33'
wnames=TRIM(wnames)//' '//'g33'

primnames= TRIM(primnames)//' '//'beta1'
wnames=TRIM(wnames)//' '//'beta1'

primnames= TRIM(primnames)//' '//'beta2'
wnames=TRIM(wnames)//' '//'beta2'

primnames= TRIM(primnames)//' '//'beta3'
wnames=TRIM(wnames)//' '//'beta3'

end subroutine specialvarnames_output
!=============================================================================
subroutine fixp_usr(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
   ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)
use mod_amrvacdef

integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision, intent(inout)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw)
double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndim)
! .. local ..
double precision, parameter        :: rhofloor = 1.0d-10
double precision, parameter        :: pfloor   = 1.0d-12
!----------------------------------------------------------------------------

where (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_) .lt. rhofloor)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_) = rhofloor
   
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,s_) = w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3,pp_) * w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)**(-eqpar(gamma_))
   
end where

where (w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,pp_) .lt. pfloor)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,pp_) = pfloor
   
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,s_) = w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3,pp_) * w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)**(-eqpar(gamma_))
   
end where

end subroutine fixp_usr
!==========================================================================================
subroutine correctaux_usr(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,patchierror,subname)

use mod_amrvacdef

integer, intent(in)            :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
integer, intent(inout)         :: patchierror(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
   ixGlo3:ixGhi3)
character(len=*), intent(in)   :: subname
double precision, intent(inout):: w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw)
double precision, intent(in)   :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndim)
! .. local ..
logical                        :: patchw(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
   ixGlo3:ixGhi3)
double precision, parameter        :: rhofloor = 1.0d-10
double precision, parameter        :: pfloor   = 1.0d-12
!-----------------------------------------------------------------------------

patchw(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3) = .true.

where (patchierror(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)/=0)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_) = rhofloor
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,pp_)  = pfloor
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,v1_)  = zero 
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,v2_)  = zero 
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,v3_)  = zero 
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,lfac_)     = 1.0d0
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,xi_)       &
      = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_) + &
      govergminone *  w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,pp_)
   patchierror(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = 0
   
   ! re-calculate entropy:
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,s_) = w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3,pp_) * w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)**(-eqpar(gamma_))
   
   patchw(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)      = .false.
end where

call conserven(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
   ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,patchw)

end subroutine correctaux_usr
!=============================================================================
subroutine flag_grid_usr(qt,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,flag)

use mod_amrvacdef

integer, intent(in)             :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
   ixGmax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
integer, intent(inout)          :: flag
double precision, intent(in)    :: qt
double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
   ixGmin3:ixGmax3,1:nw)
double precision, intent(in)    :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
   ixGmin3:ixGmax3,1:ndim)

! flag=-1 : Treat all cells active, omit deactivation (onentry, default)
! flag=0  : Treat as normal domain
! flag=1  : Treat as passive, but reduce by safety belt
! flag=2  : Always treat as passive

!-----------------------------------------------------------------------------
      
end subroutine flag_grid_usr
!=============================================================================
subroutine printlog_special

use mod_amrvacdef
!-----------------------------------------------------------------------------

call mpistop("special log file undefined")

end subroutine printlog_special
!=============================================================================
subroutine specialbound_usr(qt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iB,s)

! special boundary types, user defined
! user must assign conservative variables in bounderies

use mod_amrvacdef

integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
    ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, iB
double precision, intent(in) :: qt
type(state), intent(inout)   :: s
! .. local ..
integer                                    :: ix1,ix2,ix3, ix
integer                                    :: ixIsmin1,ixIsmin2,ixIsmin3,&
   ixIsmax1,ixIsmax2,ixIsmax3, ixIcmin1,ixIcmin2,ixIcmin3,ixIcmax1,ixIcmax2,&
   ixIcmax3
double precision                           :: wsmod(s%ws%ixGmin1:s%ws%ixGmax1,&
   s%ws%ixGmin2:s%ws%ixGmax2,s%ws%ixGmin3:s%ws%ixGmax3,1:nws)
double precision, dimension(:,:,:,:), allocatable :: xext
double precision                           :: dx1,dx2,dx3, rXmin1,rXmin2,&
   rXmin3
double precision                           :: alpha, p0
double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3)         :: r
!----------------------------------------------------------------------------
associate(x=>s%x%x,w=>s%w%w)
  
end associate
end subroutine specialbound_usr
!=============================================================================
subroutine specialrefine_grid(igrid,level,ixImin1,ixImin2,ixImin3,ixImax1,&
   ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,w,x,&
   refine,coarsen)

! Enforce additional refinement or coarsening
! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.

! you must set consistent values for integers refine/coarsen:

! refine = -1 enforce to not refine
! refine =  0 doesn't enforce anything
! refine =  1 enforce refinement

! coarsen = -1 enforce to not coarsen
! coarsen =  0 doesn't enforce anything
! coarsen =  1 enforce coarsen

use mod_amrvacdef

integer, intent(in) :: igrid, level, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision, intent(in) :: qt, w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:ndim)
integer, intent(inout) :: refine, coarsen
!-----------------------------------------------------------------------------

end subroutine specialrefine_grid
!=============================================================================
subroutine specialsource(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iwmin,iwmax,qtC,wCT,qt,w,x)

! Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
! iw=iwmin...iwmax.  wCT is at time qCT

use mod_amrvacdef

integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
    ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, iwmin,iwmax
double precision, intent(in) :: qdt, qtC, qt, x(ixImin1:ixImax1,&
   ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
double precision, intent(inout) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw), w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:nw)

! integer :: iw
! double precision :: s(ixG^T)
!-----------------------------------------------------------------------------

! do iw= iw^LIM
!    select case(iw)
!    case(m1_)
!       ! The source is based on the time centered wCT
!       call getmyforce(wCT,ixO^L,s)
!       w(ixO^S,m1_)=w(ixO^S,m1_) + qdt*s(ixO^S)
!    case(e_)
!       call getmyheating(wCT,ixO^L,s)
!       w(ixO^S,e_) =w(ixO^S,e_)  + qdt*s(ixO^S)
!    end select
! end do

end subroutine specialsource
!=============================================================================
subroutine getdt_special(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,dtnew,dx1,dx2,dx3,x)

! Limit "dt" further if necessary, e.g. due to the special source terms.
! The getdt_courant (CFL condition) and the getdt subroutine in the AMRVACPHYS
! module have already been called.

use mod_amrvacdef

integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
    ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision, intent(in) :: dx1,dx2,dx3, x(ixImin1:ixImax1,&
   ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw), dtnew
!-----------------------------------------------------------------------------

dtnew=bigdouble

end subroutine getdt_special
!=============================================================================
subroutine specialeta(w,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,idirmin,x,current,eta)

! Set the "eta" array for resistive MHD based on w or the
! "current" variable which has components between idirmin and 3.

use mod_amrvacdef

integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
    ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, idirmin
double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,nw), x(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:ndim)

double precision :: current(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
   7-2*ndir:3), eta(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
!-----------------------------------------------------------------------------

!  eta(ix^S)=...

call mpistop("con2prim can only handle constant and uniform resistivity at &
   the moment")

end subroutine specialeta
!=============================================================================
subroutine specialvarforerrest(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iflag,w,var)

! this is the place to compute a local auxiliary variable to be used
! as refinement criterion for the Lohner error estimator only
!  -->it is then requiring and iflag>nw
! note that ixO=ixI=ixG, hence the term local (gradients need special attention!)

use mod_amrvacdef

integer, intent(in)          :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iflag
double precision, intent(in) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw)
double precision, intent(out):: var(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3)
!-----------------------------------------------------------------------------

if (iflag >nw)call mpistop&
   (' iflag> nw, make change in parfile or in user file')

var(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3) = zero 

end subroutine specialvarforerrest
!=============================================================================
subroutine specialset_B0(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,wB0)

! Here one can add a steady (time-independent) potential background field

use mod_amrvacdef

integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision, intent(in)  :: x(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ixGlo3:ixGhi3,&
   1:ndim)
double precision, intent(inout) :: wB0(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndir)
!-----------------------------------------------------------------------------
call mpistop(' abs(Busr)> 0, make change in parfile or in user file')

wB0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:ndir)&
   =wB0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1:ndir)
!!wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)+user defined steady potential field

end subroutine specialset_B0
!=============================================================================

subroutine bc_int(level,qt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x)

! internal boundary, user defined
!
! This subroutine can be used to artificially overwrite ALL conservative 
! variables in a user-selected region of the mesh, and thereby act as
! an internal boundary region. It is called just before external (ghost cell)
! boundary regions will be set by the BC selection. Here, you could e.g. 
! want to introduce an extra variable (nwextra, to be distinguished from nwaux)
! which can be used to identify the internal boundary region location.
! Its effect should always be local as it acts on the mesh.
!

use mod_amrvacdef

integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
   ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,level
double precision, intent(in) :: qt
double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw)
double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndim)

! .. local ..
!logical :: patchw(ixG^T)
!----------------------------------------------------------------------------

call mpistop("bc_int not defined")

! just to give an example for relativistic MHD
!  -----------------------------------------
!patchw(ixO^S)=.true.
!where (({^D&x(ixO^S,^D)**2+})<half**2.0d0) 
!    patchw(ixO^S) = .false.
!  ^C&w(ixO^S,v^C_)=zero;
!  ^C&w(ixO^S,b^C_)=zero;
!    w(ixO^S,b3_) = one
!    w(ixO^S,v1_) = 0.99
!    w(ixO^S,rho_) = 1.d0
!    w(ixO^S,pp_)  = 2.0d0
!    w(ixO^S,lfac_)=one/dsqrt(one-({^C&w(ixO^S,v^C_)**2.0d0+}))
!end where
!!if (useprimitiveRel) then
!!  where (({^D&x(ixO^S,^D)**2+})<half**2.0d0) 
!!  {^C&w(ixO^S,u^C_)=w(ixO^S,lfac_)*w(ixO^S,v^C_);\}
!!  end where
!!endif
!call conserve(ixI^L,ixO^L,w,x,patchw)

end subroutine bc_int
!=============================================================================
subroutine process_grid_usr(igrid,level,ixImin1,ixImin2,ixImin3,ixImax1,&
   ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,qt,w,x)

! this subroutine is ONLY to be used for computing auxiliary variables
! which happen to be non-local (like div v), and are in no way used for
! flux computations. As auxiliaries, they are also not advanced

use mod_amrvacdef

integer, intent(in):: igrid,level,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision, intent(in):: qt,x(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndim)
double precision, intent(inout):: w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw)
!-----------------------------------------------------------------------------

end subroutine process_grid_usr
!=============================================================================
! amrvacusr.t
!=============================================================================
