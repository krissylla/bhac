!=============================================================================
! amrvacusr.t
!=============================================================================
!  INCLUDE:amrvacnul/speciallog.t
!  INCLUDE:amrvacnul/specialbound.t
!  INCLUDE:amrvacnul/specialsource.t
  INCLUDE:amrvacnul/specialimpl.t
!  INCLUDE:amrvacnul/usrflags.t
  INCLUDE:amrvacnul/correctaux_usr.t
  !=============================================================================
  subroutine initglobaldata_usr

    use mod_amrvacdef
    use mod_oneblock
    !-----------------------------------------------------------------------------

    eqpar(gamma_)  = 5.0d0/3.0d0 ! Adiabatic index
    
    call read_oneblock("it_53760_truncated.blk")

!    print *, woneblock(:,:,:,rho_)
!    print *, xoneblock(:,:,:,:)
  
!    call mpistop('stopping after reading')

  end subroutine initglobaldata_usr
  !=============================================================================
  subroutine initonegrid_usr(ixI^L,ixO^L,s)

    ! initialize one grid within ixO^L

    use mod_amrvacdef
    use mod_oneblock

    integer, intent(in) :: ixI^L, ixO^L
    type(state)         :: s
    ! .. local ..
    integer             :: ix^D
    
    double precision                         :: xC(ixI^S,1:ndim)
    integer                                  :: ixC^L, ixCp^L, idim, idir    
!-----------------------------------------------------------------------------
    associate(x=>s%x%x,w=>s%w%w{#IFDEF STAGGERED ,ws=>s%ws%w})

       w = 0.0d0
      
       {#IFDEF PSI
       w(ixO^S,psi_) = 0.0d0
       }

{ixCmax^D=ixOmax^D;}
{ixCmin^D=ixOmin^D-1;} ! Extend range by one    
  
!Cris attempt
!{do ix^D=ixOmin^D, ixOmax^D \}
!    {#IFNDEF STAGGERED
!    call interpolate_oneblock( x(ix^D,:) , b1_, w(ix^D, b1_) )
!    call interpolate_oneblock( x(ix^D,:) , b2_, w(ix^D, b2_) )
!    call interpolate_oneblock( x(ix^D,:) , b3_, w(ix^D, b3_) )}
!    {#IFDEF STAGGERED
!    call interpolate_oneblock( x(ix^D,:) , b1_, ws(ix^D, bs1_) )
!    call interpolate_oneblock( x(ix^D,:) , b2_, ws(ix^D, bs2_) )
!    call interpolate_oneblock( x(ix^D,:) , b3_, ws(ix^D, bs3_) )}
!    !print *, x(ix^D,:), ws(ix^D, b1_)
!{enddo\}

do idim=1,ndim
          ! Get edge coordinates
          do idir=1,ndim
            if (idim.eq.idir) then
              {ixCp^L=ixC^L+kr(idir,^D);}
              xC(ixC^S,idir) = half * (x(ixCp^S,idir) + x(ixC^S,idir))
            else
              xC(ixC^S,idir)=x(ixC^S,idir)
            end if
          end do

  select case(idim)
  case(1)
  {do ix^D=ixCmin^D, ixCmax^D \}
    call interpolate_oneblock( xC(ix^D,:) , b1_, ws(ix^D, bs1_) )
  {enddo\}

  case(2)
  {do ix^D=ixCmin^D, ixCmax^D \}
    call interpolate_oneblock( xC(ix^D,:) , b2_, ws(ix^D, bs2_) )
  {enddo\}

  case(3)
  {do ix^D=ixCmin^D, ixCmax^D \}
    call interpolate_oneblock( xC(ix^D,:) , b3_, ws(ix^D, bs3_) )
  {enddo\}

  end select
end do


call faces2centers(ixO^L,s)

{do ix^D=ixOmin^D, ixOmax^D \}
    !call interpolate_oneblock( x(ix^D,:) , b3_, w(ix^D, b3_) ) 
    call interpolate_oneblock( x(ix^D,:) , rho_, w(ix^D, rho_) )
    call interpolate_oneblock( x(ix^D,:) , pp_, w(ix^D, pp_) )

    call interpolate_oneblock( x(ix^D,:) , u1_, w(ix^D, u1_) )
    call interpolate_oneblock( x(ix^D,:) , u2_, w(ix^D, u2_) )
    call interpolate_oneblock( x(ix^D,:) , u3_, w(ix^D, u3_) )
    
    

{enddo\}

       
      {#IFDEF ENTROPY
      w(ixO^S,s_) = w(ixO^S,pp_) * w(ixO^S,rho_)**(-eqpar(gamma_))
      }

    call conserve(ixI^L,ixO^L,w,x,patchfalse)
 
  end associate
end subroutine initonegrid_usr
!=============================================================================
subroutine initvecpot_usr(ixI^L, ixC^L, xC, A, idir)

  ! initialize the vectorpotential on the corners
  ! used by b_from_vectorpotential()


  use mod_amrvacdef

  integer, intent(in)                :: ixI^L, ixC^L, idir
  double precision, intent(in)       :: xC(ixI^S,1:ndim)
  double precision, intent(out)      :: A(ixI^S)
  ! .. local ..
  !-----------------------------------------------------------------------------

end subroutine initvecpot_usr
!=============================================================================
subroutine specialvar_output(ixI^L,ixO^L,nwmax,w,s,normconv)

! this subroutine can be used in convert, to add auxiliary variables to the
! converted output file, for further analysis using tecplot, paraview, ....
! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
!
! the array normconv can be filled in the (nw+1:nw+nwauxio) range with 
! corresponding normalization values (default value 1)

use mod_amrvacdef

integer, intent(in)                :: ixI^L,ixO^L,nwmax
double precision                   :: w(ixI^S,1:nwmax)
type(state)                        :: s
double precision                   :: normconv(0:nwmax)
!-----------------------------------------------------------------------------
associate(x=>s%x%x{#IFDEF STAGGERED ,ws=>s%ws%w})

{#IFDEF STAGGERED
call div_staggered(ixO^L,s,w(ixO^S,nw+1))
}{#IFNDEF STAGGERED
! Reduce output array size, +1 was added for eventual pointdata output
call get_divb(ixI^L,ixO^L^LSUB1,w(ixI^S,1:nw),w(ixI^S,nw+1))
}

end associate
end subroutine specialvar_output
!=============================================================================
subroutine fixp_usr(ixI^L,ixO^L,w,x)
use mod_amrvacdef

integer, intent(in)                :: ixI^L,ixO^L
double precision, intent(inout)    :: w(ixI^S,1:nw)
double precision, intent(in)       :: x(ixI^S,1:ndim)
! .. local ..
double precision, parameter        :: rhofloor=1.0d-15
!----------------------------------------------------------------------------

where (w(ixO^S,rho_) .lt. rhofloor)
   w(ixO^S,rho_) = rhofloor
   {#IFDEF ENTROPY
   w(ixO^S,s_) = w(ixO^S,pp_) * w(ixO^S,rho_)**(-eqpar(gamma_))
   }
end where

end subroutine fixp_usr
!=============================================================================
subroutine flag_grid_usr(qt,ixG^L,ixO^L,w,x,flag)

use mod_amrvacdef

integer, intent(in)             :: ixG^L, ixO^L
integer, intent(inout)          :: flag
double precision, intent(in)    :: qt
double precision, intent(inout) :: w(ixG^S,1:nw)
double precision, intent(in)    :: x(ixG^S,1:ndim)

! flag=-1 : Treat all cells active, omit deactivation (onentry, default)
! flag=0  : Treat as normal domain
! flag=1  : Treat as passive, but reduce by safety belt
! flag=2  : Always treat as passive

!-----------------------------------------------------------------------------
      
end subroutine flag_grid_usr
!=============================================================================
subroutine specialvarnames_output

! newly added variables to be concatenated with the primnames/wnames string

use mod_amrvacdef
!-----------------------------------------------------------------------------

primnames= TRIM(primnames)//' '//'divB'
wnames=TRIM(wnames)//' '//'divB'

end subroutine specialvarnames_output
!=============================================================================
subroutine printlog_special

use mod_amrvacdef
!-----------------------------------------------------------------------------

call mpistop("special log file undefined")

end subroutine printlog_special
!=============================================================================
subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,s)

! special boundary types, user defined
! user must assign conservative variables in bounderies

use mod_amrvacdef

integer, intent(in) :: ixI^L, ixO^L, iB
double precision, intent(in) :: qt
type(state), intent(inout)   :: s
! .. local ..
integer                                    :: ix^D, ix
integer                                    :: ixIs^L, ixIc^L
double precision                           :: wsmod(s%ws%ixG^S,1:nws)
double precision, dimension(:,:,:,:), allocatable :: xext
double precision                           :: dx^D, rXmin^D
double precision                           :: alpha, p0
double precision, dimension(ixI^S)         :: r
!----------------------------------------------------------------------------
associate(x=>s%x%x,w=>s%w%w{#IFDEF STAGGERED ,ws=>s%ws%w})
  
end associate
end subroutine specialbound_usr
!=============================================================================
subroutine specialrefine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)

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

integer, intent(in) :: igrid, level, ixI^L, ixO^L
double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
integer, intent(inout) :: refine, coarsen
!-----------------------------------------------------------------------------

end subroutine specialrefine_grid
!=============================================================================
subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

! Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
! iw=iwmin...iwmax.  wCT is at time qCT

use mod_amrvacdef

integer, intent(in) :: ixI^L, ixO^L, iw^LIM
double precision, intent(in) :: qdt, qtC, qt, x(ixI^S,1:ndim)
double precision, intent(inout) :: wCT(ixI^S,1:nw), w(ixI^S,1:nw)

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
subroutine getdt_special(w,ixI^L,ixO^L,dtnew,dx^D,x)

! Limit "dt" further if necessary, e.g. due to the special source terms.
! The getdt_courant (CFL condition) and the getdt subroutine in the AMRVACPHYS
! module have already been called.

use mod_amrvacdef

integer, intent(in) :: ixI^L, ixO^L
double precision, intent(in) :: dx^D, x(ixI^S,1:ndim)
double precision, intent(inout) :: w(ixI^S,1:nw), dtnew
!-----------------------------------------------------------------------------

dtnew=bigdouble

end subroutine getdt_special
!=============================================================================
subroutine specialeta(w,ixI^L,ixO^L,idirmin,x,current,eta)

! Set the "eta" array for resistive MHD based on w or the
! "current" variable which has components between idirmin and 3.

use mod_amrvacdef

integer, intent(in) :: ixI^L, ixO^L, idirmin
double precision, intent(in) :: w(ixI^S,nw), x(ixI^S,1:ndim)

double precision :: current(ixG^T,7-2*ndir:3), eta(ixG^T)
!-----------------------------------------------------------------------------

!  eta(ix^S)=...

call mpistop("con2prim can only handle constant and uniform resistivity at the moment")

end subroutine specialeta
!=============================================================================
subroutine specialvarforerrest(ixI^L,ixO^L,iflag,w,var)

! this is the place to compute a local auxiliary variable to be used
! as refinement criterion for the Lohner error estimator only
!  -->it is then requiring and iflag>nw
! note that ixO=ixI=ixG, hence the term local (gradients need special attention!)

use mod_amrvacdef

integer, intent(in)          :: ixI^L,ixO^L,iflag
double precision, intent(in) :: w(ixI^S,1:nw)
double precision, intent(out):: var(ixG^T)
!-----------------------------------------------------------------------------

if (iflag >nw)call mpistop(' iflag> nw, make change in parfile or in user file')

var(ixI^S) = zero 

end subroutine specialvarforerrest
!=============================================================================
subroutine specialset_B0(ixI^L,ixO^L,x,wB0)

! Here one can add a steady (time-independent) potential background field

use mod_amrvacdef

integer, intent(in)           :: ixI^L,ixO^L
double precision, intent(in)  :: x(ixG^T,1:ndim)
double precision, intent(inout) :: wB0(ixI^S,1:ndir)
!-----------------------------------------------------------------------------
call mpistop(' abs(Busr)> 0, make change in parfile or in user file')

wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)
!!wB0(ixO^S,1:ndir)=wB0(ixO^S,1:ndir)+user defined steady potential field

end subroutine specialset_B0
!=============================================================================

subroutine bc_int(level,qt,ixI^L,ixO^L,w,x)

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
use mod_oneblock

integer, intent(in) :: ixI^L,ixO^L,level
double precision, intent(in) :: qt
double precision, intent(inout) :: w(ixI^S,1:nw)
double precision, intent(in) :: x(ixI^S,1:ndim)
integer             :: ix^D

! .. local ..
logical :: patchw(ixG^T)
!----------------------------------------------------------------------------

!call mpistop("bc_int not defined")

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

!Cris attempt

patchw = .true.
{do ix^D=ixOmin^D, ixOmax^D \}
	if (x(ix^D,3)<22.0d0) then
	patchw(ix^D) = .false.
    
    call interpolate_oneblock( x(ix^D,:) , rho_, w(ix^D, rho_) )
    call interpolate_oneblock( x(ix^D,:) , pp_, w(ix^D, pp_) )

    call interpolate_oneblock( x(ix^D,:) , u1_, w(ix^D, u1_) )
    call interpolate_oneblock( x(ix^D,:) , u2_, w(ix^D, u2_) )
    call interpolate_oneblock( x(ix^D,:) , u3_, w(ix^D, u3_) )
	end if
{enddo\}

{#IFDEF ENTROPY
{do ix^D=ixOmin^D, ixOmax^D \}
		if (x(ix^D,3)<22.0d0) then
      w(ix^D,s_) = w(ix^D,pp_) * w(ix^D,rho_)**(-eqpar(gamma_))      
      end if
{enddo\}}

 call conserve(ixI^L,ixO^L,w,x,patchw)


end subroutine bc_int
!=============================================================================
subroutine process_grid_usr(igrid,level,ixI^L,ixO^L,qt,w,x)

! this subroutine is ONLY to be used for computing auxiliary variables
! which happen to be non-local (like div v), and are in no way used for
! flux computations. As auxiliaries, they are also not advanced

use mod_amrvacdef

integer, intent(in):: igrid,level,ixI^L,ixO^L
double precision, intent(in):: qt,x(ixI^S,1:ndim)
double precision, intent(inout):: w(ixI^S,1:nw)
!-----------------------------------------------------------------------------

end subroutine process_grid_usr
!=============================================================================
! amrvacusr.t
!=============================================================================
