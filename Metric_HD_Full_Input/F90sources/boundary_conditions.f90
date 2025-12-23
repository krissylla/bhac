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
subroutine bc_phys(iside,idims,time,s,ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,&
   ixBmax3)

use mod_amrvacdef

integer, intent(in)          :: iside, idims, ixBmin1,ixBmin2,ixBmin3,ixBmax1,&
   ixBmax2,ixBmax3
double precision, intent(in) :: time
type(state), intent(inout)   :: s

integer :: iw, iB, ix1,ix2,ix3, ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3, ixDmin1,ixDmin2,ixDmin3,ixDmax1,ixDmax2,ixDmax3, ixIpmin1,&
   ixIpmin2,ixIpmin3,ixIpmax1,ixIpmax2,ixIpmax3, ixGmin1,ixGmin2,ixGmin3,&
   ixGmax1,ixGmax2,ixGmax3

double precision, dimension(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
   1:nw)    :: wsave
integer, parameter :: dixPolefix=2, dixHotaka=1
logical            :: is_coarse(ndim,2)
!-----------------------------------------------------------------------------
associate(x=>s%x%x,w=>s%w%w)
ixGmin1=s%w%ixGmin1;ixGmin2=s%w%ixGmin2;ixGmin3=s%w%ixGmin3
ixGmax1=s%w%ixGmax1;ixGmax2=s%w%ixGmax2;ixGmax3=s%w%ixGmax3;

! Test if we have a coarse neighbor:
call is_neighbor_coarse(s,is_coarse)

select case (idims)
case (1)
if (iside==2) then

   
   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! maximal boundary
   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      iB=ismax1
      ixImin1=ixBmax1+1-dixB;ixImin2=ixBmin2;ixImin3=ixBmin3;
      ixImax1=ixBmax1;ixImax2=ixBmax2;ixImax3=ixBmax3;

      !!!!!!! Make tangential ranges greater for staggered components
      !!!!!!! Add fill the normal component

      ! Apply on primitive variables?
      if (primitiveB(iside,idims)) then
         ixDmin1=ixBmin1+dixB;ixDmin2=ixBmin2;ixDmin3=ixBmin3;
         ixDmax1=ixBmax1-dixB;ixDmax2=ixBmax2;ixDmax3=ixBmax3;
         wsave(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) &
            = w(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) !make a backup of the state in domain since p2c/c2p might lead to slight changes in w
         call primitive(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
            ixDmin1,ixDmin2,ixDmin3,ixDmax1,ixDmax2,ixDmax3,w,x)
      end if
      
      ! cont/symm/asymm types
      do iw=1,nwflux+nwaux

         select case (typeB(iw,iB))
         case ("symm")
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               = w(ixImin1-1:ixImin1-dixB:-1,ixImin2:ixImax2,ixImin3:ixImax3,&
               iw)
         case ("asymm")
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               =-w(ixImin1-1:ixImin1-dixB:-1,ixImin2:ixImax2,ixImin3:ixImax3,&
               iw)
         case ("cont")
            do ix1=ixImin1,ixImax1
               w(ix1,ixImin2:ixImax2,ixImin3:ixImax3,iw) = w(ixImin1-1,&
                  ixImin2:ixImax2,ixImin3:ixImax3,iw)
            end do
         case("noinflow")
            if (iw==1+1)then
              do ix1=ixImin1,ixImax1
                  w(ix1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
                     = max(w(ixImin1-1,ixImin2:ixImax2,ixImin3:ixImax3,iw),&
                     zero)
              end do
            else
              do ix1=ixImin1,ixImax1
                  w(ix1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
                     = w(ixImin1-1,ixImin2:ixImax2,ixImin3:ixImax3,iw)
              end do
            end if
         case("limitinflow")
            if (iw==1+1)then
              do ix1=ixImin1,ixImax1
                  w(ix1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
                     = max(w(ixImin1-1,ixImin2:ixImax2,ixImin3:ixImax3,iw), &
                                           w(ixImin1-1,ixImin2:ixImax2,&
                                              ixImin3:ixImax3,iw)*ratebdflux)
              end do
            else
              do ix1=ixImin1,ixImax1
                  w(ix1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
                     = w(ixImin1-1,ixImin2:ixImax2,ixImin3:ixImax3,iw)
              end do
           end if
        case ("polefix")
           ! First overwrite cells in domain
           ixIpmin1=ixImin1-dixPolefix;ixIpmin2=ixImin2;ixIpmin3=ixImin3;
           ixIpmax1=ixImin1-1;ixIpmax2=ixImax2;ixIpmax3=ixImax3;
           do ix1=ixIpmin1,ixIpmax1
              w(ix1,ixIpmin2:ixIpmax2,ixIpmin3:ixIpmax3,iw) &
                 = w(ixIpmin1-1,ixIpmin2:ixIpmax2,ixIpmin3:ixIpmax3,iw)
           end do
           ! Now apply symm boundary
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               = w(ixImin1-1:ixImin1-dixB:-1,ixImin2:ixImax2,ixImin3:ixImax3,&
               iw)
        case ("apolefix")
           ! First interpolate cells in domain down to zero at pole
           ixIpmin1=ixImin1-dixPolefix;ixIpmin2=ixImin2;ixIpmin3=ixImin3;
           ixIpmax1=ixImax1-1;ixIpmax2=ixImax2;ixIpmax3=ixImax3;
           do ix1=ixIpmin1,ixIpmax1
              w(ix1,ixIpmin2:ixIpmax2,ixIpmin3:ixIpmax3,iw) &
                 = w(ixIpmin1-1,ixIpmin2:ixIpmax2,ixIpmin3:ixIpmax3,iw) &
                   * (xprobmax1-x(ix1,ixIpmin2:ixIpmax2,ixIpmin3:ixIpmax3,&
                      1)) &
                   / (xprobmax1-x(ixIpmin1-1,ixIpmin2:ixIpmax2,&
                      ixIpmin3:ixIpmax3,1))
           end do
        case ("hotaka")
           ! First overwrite cells in domain with zero
           ixIpmin1=ixImin1-dixHotaka;ixIpmin2=ixImin2;ixIpmin3=ixImin3;
           ixIpmax1=ixImin1-1;ixIpmax2=ixImax2;ixIpmax3=ixImax3;
           do ix1=ixIpmin1,ixIpmax1
              w(ix1,ixIpmin2:ixIpmax2,ixIpmin3:ixIpmax3,iw) = zero
           end do
           ! Now apply symm boundary
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               = w(ixImin1-1:ixImin1-dixB:-1,ixImin2:ixImax2,ixImin3:ixImax3,&
               iw)
        case ("ahotaka")
            ! First overwrite cells in domain with zero
           ixIpmin1=ixImin1-dixHotaka;ixIpmin2=ixImin2;ixIpmin3=ixImin3;
           ixIpmax1=ixImax1-1;ixIpmax2=ixImax2;ixIpmax3=ixImax3;
           do ix1=ixIpmin1,ixIpmax1
              w(ix1,ixIpmin2:ixIpmax2,ixIpmin3:ixIpmax3,iw) = zero
           end do
           ! Now apply asymm boundary
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               =-w(ixImin1-1:ixImin1-dixB:-1,ixImin2:ixImax2,ixImin3:ixImax3,&
               iw)
         case ("special")
            ! skip it here, do AFTER all normal type boundaries are set
         case ("aperiodic")
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               = - w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw)
         case ("periodic")
!            call mpistop("periodic bc info should come from neighbors")
         case default
            write (unitterm,*) "Undefined boundarytype ",typeB(iw,iB), &
               "for variable iw=",iw," and side iB=",iB
         end select


      end do


      ! Apply on primitive variables?
      if (primitiveB(iside,idims)) then
         w(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) &
            = wsave(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) !restore from backup of the state in domain since p2c/c2p might lead to slight changes in w
         call conserve(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
            ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,w,x,patchfalse)
      end if
      
   else
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! minimal boundary
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      iB=ismin1
      ixImin1=ixBmin1;ixImin2=ixBmin2;ixImin3=ixBmin3;
      ixImax1=ixBmin1-1+dixB;ixImax2=ixBmax2;ixImax3=ixBmax3;

      ! Apply on primitive variables?
      if (primitiveB(iside,idims)) then         
         ixDmin1=ixBmin1+dixB;ixDmin2=ixBmin2;ixDmin3=ixBmin3;
         ixDmax1=ixBmax1-dixB;ixDmax2=ixBmax2;ixDmax3=ixBmax3;
         wsave(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) &
            = w(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) !make a backup of the state in domain since p2c/c2p might lead to slight changes in w
         call primitive(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
            ixDmin1,ixDmin2,ixDmin3,ixDmax1,ixDmax2,ixDmax3,w,x)
      end if
 
      ! cont/symm/asymm types
      do iw=1,nwflux+nwaux

   select case (typeB(iw,iB))
   case ("symm")
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         = w(ixImax1+dixB:ixImax1+1:-1,ixImin2:ixImax2,ixImin3:ixImax3,iw)
   case ("asymm")
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         =-w(ixImax1+dixB:ixImax1+1:-1,ixImin2:ixImax2,ixImin3:ixImax3,iw)
   case ("cont")
      do ix1=ixImin1,ixImax1
         w(ix1,ixImin2:ixImax2,ixImin3:ixImax3,iw) = w(ixImax1+1,&
            ixImin2:ixImax2,ixImin3:ixImax3,iw)
      end do
   case("noinflow")
      if (iw==1+1)then
         do ix1=ixImin1,ixImax1
            w(ix1,ixImin2:ixImax2,ixImin3:ixImax3,iw) = min(w(ixImax1+1,&
               ixImin2:ixImax2,ixImin3:ixImax3,iw),zero)
         end do
      else
         do ix1=ixImin1,ixImax1
            w(ix1,ixImin2:ixImax2,ixImin3:ixImax3,iw) = w(ixImax1+1,&
               ixImin2:ixImax2,ixImin3:ixImax3,iw)
         end do
      end if
   case("limitinflow")
      if (iw==1+1)then
         do ix1=ixImin1,ixImax1
            w(ix1,ixImin2:ixImax2,ixImin3:ixImax3,iw) = min(w(ixImax1+1,&
               ixImin2:ixImax2,ixImin3:ixImax3,iw), &
                 w(ixImax1+1,ixImin2:ixImax2,ixImin3:ixImax3,iw)*ratebdflux)
         end do
      else
         do ix1=ixImin1,ixImax1
            w(ix1,ixImin2:ixImax2,ixImin3:ixImax3,iw) = w(ixImax1+1,&
               ixImin2:ixImax2,ixImin3:ixImax3,iw)
         end do
      end if
   case ("polefix")
      ! First overwrite cells in domain
      ixIpmin1=ixImax1+1;ixIpmin2=ixImin2;ixIpmin3=ixImin3;
      ixIpmax1=ixImax1+dixPolefix;ixIpmax2=ixImax2;ixIpmax3=ixImax3;
      do ix1=ixIpmin1,ixIpmax1
         w(ix1,ixIpmin2:ixIpmax2,ixIpmin3:ixIpmax3,iw) = w(ixIpmax1+1,&
            ixIpmin2:ixIpmax2,ixIpmin3:ixIpmax3,iw)
      end do
      ! Now apply symm boundary
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         = w(ixImax1+dixB:ixImax1+1:-1,ixImin2:ixImax2,ixImin3:ixImax3,iw)
   case ("apolefix")
      ! First interpolate cells in domain down to zero at pole
      ixIpmin1=ixImax1+1;ixIpmin2=ixImin2;ixIpmin3=ixImin3;
      ixIpmax1=ixImax1+dixPolefix;ixIpmax2=ixImax2;ixIpmax3=ixImax3;
      do ix1=ixIpmin1,ixIpmax1
         w(ix1,ixIpmin2:ixIpmax2,ixIpmin3:ixIpmax3,iw) = w(ixIpmax1+1,&
            ixIpmin2:ixIpmax2,ixIpmin3:ixIpmax3,iw) &
              * (x(ix1,ixIpmin2:ixIpmax2,ixIpmin3:ixIpmax3,1))&
                 /x(ixIpmax1+1,ixIpmin2:ixIpmax2,ixIpmin3:ixIpmax3,1)
      end do
      ! Now apply asymm boundary
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         =-w(ixImax1+dixB:ixImax1+1:-1,ixImin2:ixImax2,ixImin3:ixImax3,iw)
   case ("hotaka")
      ! First overwrite cells in domain with zero
      ixIpmin1=ixImax1+1;ixIpmin2=ixImin2;ixIpmin3=ixImin3;
      ixIpmax1=ixImax1+dixHotaka;ixIpmax2=ixImax2;ixIpmax3=ixImax3;
      do ix1=ixIpmin1,ixIpmax1
         w(ix1,ixIpmin2:ixIpmax2,ixIpmin3:ixIpmax3,iw) = zero
      end do
      ! Now apply symm boundary
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         = w(ixImax1+dixB:ixImax1+1:-1,ixImin2:ixImax2,ixImin3:ixImax3,iw)
   case ("ahotaka")
      ! First overwrite cells in domain with zero
      ixIpmin1=ixImax1+1;ixIpmin2=ixImin2;ixIpmin3=ixImin3;
      ixIpmax1=ixImax1+dixHotaka;ixIpmax2=ixImax2;ixIpmax3=ixImax3;
      do ix1=ixIpmin1,ixIpmax1
         w(ix1,ixIpmin2:ixIpmax2,ixIpmin3:ixIpmax3,iw) = zero
      end do
      ! Now apply asymm boundary
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         =-w(ixImax1+dixB:ixImax1+1:-1,ixImin2:ixImax2,ixImin3:ixImax3,iw)
   case ("special")
      ! skip it here, do AFTER all normal type boundaries are set
   case ("aperiodic")
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         = - w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw)
   case ("periodic")
      !            call mpistop("periodic bc info should come from neighbors")
   case default
      write (unitterm,*) "Undefined boundarytype ",typeB(iw,iB), &
           "for variable iw=",iw," and side iB=",iB
   end select


end do


      
      ! Apply on primitive variables?
      if (primitiveB(iside,idims)) then
         w(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) &
            = wsave(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) !restore from backup of the state in domain since p2c/c2p might lead to slight changes in w
         call conserve(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
            ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,w,x,patchfalse)
      end if
      
   end if 
case (2)
if (iside==2) then

   
   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! maximal boundary
   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      iB=ismax2
      ixImin1=ixBmin1;ixImin2=ixBmax2+1-dixB;ixImin3=ixBmin3;
      ixImax1=ixBmax1;ixImax2=ixBmax2;ixImax3=ixBmax3;

      !!!!!!! Make tangential ranges greater for staggered components
      !!!!!!! Add fill the normal component

      ! Apply on primitive variables?
      if (primitiveB(iside,idims)) then
         ixDmin1=ixBmin1;ixDmin2=ixBmin2+dixB;ixDmin3=ixBmin3;
         ixDmax1=ixBmax1;ixDmax2=ixBmax2-dixB;ixDmax3=ixBmax3;
         wsave(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) &
            = w(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) !make a backup of the state in domain since p2c/c2p might lead to slight changes in w
         call primitive(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
            ixDmin1,ixDmin2,ixDmin3,ixDmax1,ixDmax2,ixDmax3,w,x)
      end if
      
      ! cont/symm/asymm types
      do iw=1,nwflux+nwaux

         select case (typeB(iw,iB))
         case ("symm")
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               = w(ixImin1:ixImax1,ixImin2-1:ixImin2-dixB:-1,ixImin3:ixImax3,&
               iw)
         case ("asymm")
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               =-w(ixImin1:ixImax1,ixImin2-1:ixImin2-dixB:-1,ixImin3:ixImax3,&
               iw)
         case ("cont")
            do ix2=ixImin2,ixImax2
               w(ixImin1:ixImax1,ix2,ixImin3:ixImax3,iw) = w(ixImin1:ixImax1,&
                  ixImin2-1,ixImin3:ixImax3,iw)
            end do
         case("noinflow")
            if (iw==1+2)then
              do ix2=ixImin2,ixImax2
                  w(ixImin1:ixImax1,ix2,ixImin3:ixImax3,iw) &
                     = max(w(ixImin1:ixImax1,ixImin2-1,ixImin3:ixImax3,iw),&
                     zero)
              end do
            else
              do ix2=ixImin2,ixImax2
                  w(ixImin1:ixImax1,ix2,ixImin3:ixImax3,iw) &
                     = w(ixImin1:ixImax1,ixImin2-1,ixImin3:ixImax3,iw)
              end do
            end if
         case("limitinflow")
            if (iw==1+2)then
              do ix2=ixImin2,ixImax2
                  w(ixImin1:ixImax1,ix2,ixImin3:ixImax3,iw) &
                     = max(w(ixImin1:ixImax1,ixImin2-1,ixImin3:ixImax3,iw), &
                                           w(ixImin1:ixImax1,ixImin2-1,&
                                              ixImin3:ixImax3,iw)*ratebdflux)
              end do
            else
              do ix2=ixImin2,ixImax2
                  w(ixImin1:ixImax1,ix2,ixImin3:ixImax3,iw) &
                     = w(ixImin1:ixImax1,ixImin2-1,ixImin3:ixImax3,iw)
              end do
           end if
        case ("polefix")
           ! First overwrite cells in domain
           ixIpmin1=ixImin1;ixIpmin2=ixImin2-dixPolefix;ixIpmin3=ixImin3;
           ixIpmax1=ixImax1;ixIpmax2=ixImin2-1;ixIpmax3=ixImax3;
           do ix2=ixIpmin2,ixIpmax2
              w(ixIpmin1:ixIpmax1,ix2,ixIpmin3:ixIpmax3,iw) &
                 = w(ixIpmin1:ixIpmax1,ixIpmin2-1,ixIpmin3:ixIpmax3,iw)
           end do
           ! Now apply symm boundary
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               = w(ixImin1:ixImax1,ixImin2-1:ixImin2-dixB:-1,ixImin3:ixImax3,&
               iw)
        case ("apolefix")
           ! First interpolate cells in domain down to zero at pole
           ixIpmin1=ixImin1;ixIpmin2=ixImin2-dixPolefix;ixIpmin3=ixImin3;
           ixIpmax1=ixImax1;ixIpmax2=ixImax2-1;ixIpmax3=ixImax3;
           do ix2=ixIpmin2,ixIpmax2
              w(ixIpmin1:ixIpmax1,ix2,ixIpmin3:ixIpmax3,iw) &
                 = w(ixIpmin1:ixIpmax1,ixIpmin2-1,ixIpmin3:ixIpmax3,iw) &
                   * (xprobmax2-x(ixIpmin1:ixIpmax1,ix2,ixIpmin3:ixIpmax3,&
                      2)) &
                   / (xprobmax2-x(ixIpmin1:ixIpmax1,ixIpmin2-1,&
                      ixIpmin3:ixIpmax3,2))
           end do
        case ("hotaka")
           ! First overwrite cells in domain with zero
           ixIpmin1=ixImin1;ixIpmin2=ixImin2-dixHotaka;ixIpmin3=ixImin3;
           ixIpmax1=ixImax1;ixIpmax2=ixImin2-1;ixIpmax3=ixImax3;
           do ix2=ixIpmin2,ixIpmax2
              w(ixIpmin1:ixIpmax1,ix2,ixIpmin3:ixIpmax3,iw) = zero
           end do
           ! Now apply symm boundary
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               = w(ixImin1:ixImax1,ixImin2-1:ixImin2-dixB:-1,ixImin3:ixImax3,&
               iw)
        case ("ahotaka")
            ! First overwrite cells in domain with zero
           ixIpmin1=ixImin1;ixIpmin2=ixImin2-dixHotaka;ixIpmin3=ixImin3;
           ixIpmax1=ixImax1;ixIpmax2=ixImax2-1;ixIpmax3=ixImax3;
           do ix2=ixIpmin2,ixIpmax2
              w(ixIpmin1:ixIpmax1,ix2,ixIpmin3:ixIpmax3,iw) = zero
           end do
           ! Now apply asymm boundary
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               =-w(ixImin1:ixImax1,ixImin2-1:ixImin2-dixB:-1,ixImin3:ixImax3,&
               iw)
         case ("special")
            ! skip it here, do AFTER all normal type boundaries are set
         case ("aperiodic")
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               = - w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw)
         case ("periodic")
!            call mpistop("periodic bc info should come from neighbors")
         case default
            write (unitterm,*) "Undefined boundarytype ",typeB(iw,iB), &
               "for variable iw=",iw," and side iB=",iB
         end select


      end do


      ! Apply on primitive variables?
      if (primitiveB(iside,idims)) then
         w(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) &
            = wsave(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) !restore from backup of the state in domain since p2c/c2p might lead to slight changes in w
         call conserve(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
            ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,w,x,patchfalse)
      end if
      
   else
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! minimal boundary
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      iB=ismin2
      ixImin1=ixBmin1;ixImin2=ixBmin2;ixImin3=ixBmin3;
      ixImax1=ixBmax1;ixImax2=ixBmin2-1+dixB;ixImax3=ixBmax3;

      ! Apply on primitive variables?
      if (primitiveB(iside,idims)) then         
         ixDmin1=ixBmin1;ixDmin2=ixBmin2+dixB;ixDmin3=ixBmin3;
         ixDmax1=ixBmax1;ixDmax2=ixBmax2-dixB;ixDmax3=ixBmax3;
         wsave(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) &
            = w(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) !make a backup of the state in domain since p2c/c2p might lead to slight changes in w
         call primitive(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
            ixDmin1,ixDmin2,ixDmin3,ixDmax1,ixDmax2,ixDmax3,w,x)
      end if
 
      ! cont/symm/asymm types
      do iw=1,nwflux+nwaux

   select case (typeB(iw,iB))
   case ("symm")
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         = w(ixImin1:ixImax1,ixImax2+dixB:ixImax2+1:-1,ixImin3:ixImax3,iw)
   case ("asymm")
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         =-w(ixImin1:ixImax1,ixImax2+dixB:ixImax2+1:-1,ixImin3:ixImax3,iw)
   case ("cont")
      do ix2=ixImin2,ixImax2
         w(ixImin1:ixImax1,ix2,ixImin3:ixImax3,iw) = w(ixImin1:ixImax1,&
            ixImax2+1,ixImin3:ixImax3,iw)
      end do
   case("noinflow")
      if (iw==1+2)then
         do ix2=ixImin2,ixImax2
            w(ixImin1:ixImax1,ix2,ixImin3:ixImax3,iw) = min(w(ixImin1:ixImax1,&
               ixImax2+1,ixImin3:ixImax3,iw),zero)
         end do
      else
         do ix2=ixImin2,ixImax2
            w(ixImin1:ixImax1,ix2,ixImin3:ixImax3,iw) = w(ixImin1:ixImax1,&
               ixImax2+1,ixImin3:ixImax3,iw)
         end do
      end if
   case("limitinflow")
      if (iw==1+2)then
         do ix2=ixImin2,ixImax2
            w(ixImin1:ixImax1,ix2,ixImin3:ixImax3,iw) = min(w(ixImin1:ixImax1,&
               ixImax2+1,ixImin3:ixImax3,iw), &
                 w(ixImin1:ixImax1,ixImax2+1,ixImin3:ixImax3,iw)*ratebdflux)
         end do
      else
         do ix2=ixImin2,ixImax2
            w(ixImin1:ixImax1,ix2,ixImin3:ixImax3,iw) = w(ixImin1:ixImax1,&
               ixImax2+1,ixImin3:ixImax3,iw)
         end do
      end if
   case ("polefix")
      ! First overwrite cells in domain
      ixIpmin1=ixImin1;ixIpmin2=ixImax2+1;ixIpmin3=ixImin3;
      ixIpmax1=ixImax1;ixIpmax2=ixImax2+dixPolefix;ixIpmax3=ixImax3;
      do ix2=ixIpmin2,ixIpmax2
         w(ixIpmin1:ixIpmax1,ix2,ixIpmin3:ixIpmax3,iw) = w(ixIpmin1:ixIpmax1,&
            ixIpmax2+1,ixIpmin3:ixIpmax3,iw)
      end do
      ! Now apply symm boundary
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         = w(ixImin1:ixImax1,ixImax2+dixB:ixImax2+1:-1,ixImin3:ixImax3,iw)
   case ("apolefix")
      ! First interpolate cells in domain down to zero at pole
      ixIpmin1=ixImin1;ixIpmin2=ixImax2+1;ixIpmin3=ixImin3;
      ixIpmax1=ixImax1;ixIpmax2=ixImax2+dixPolefix;ixIpmax3=ixImax3;
      do ix2=ixIpmin2,ixIpmax2
         w(ixIpmin1:ixIpmax1,ix2,ixIpmin3:ixIpmax3,iw) = w(ixIpmin1:ixIpmax1,&
            ixIpmax2+1,ixIpmin3:ixIpmax3,iw) &
              * (x(ixIpmin1:ixIpmax1,ix2,ixIpmin3:ixIpmax3,2))&
                 /x(ixIpmin1:ixIpmax1,ixIpmax2+1,ixIpmin3:ixIpmax3,2)
      end do
      ! Now apply asymm boundary
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         =-w(ixImin1:ixImax1,ixImax2+dixB:ixImax2+1:-1,ixImin3:ixImax3,iw)
   case ("hotaka")
      ! First overwrite cells in domain with zero
      ixIpmin1=ixImin1;ixIpmin2=ixImax2+1;ixIpmin3=ixImin3;
      ixIpmax1=ixImax1;ixIpmax2=ixImax2+dixHotaka;ixIpmax3=ixImax3;
      do ix2=ixIpmin2,ixIpmax2
         w(ixIpmin1:ixIpmax1,ix2,ixIpmin3:ixIpmax3,iw) = zero
      end do
      ! Now apply symm boundary
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         = w(ixImin1:ixImax1,ixImax2+dixB:ixImax2+1:-1,ixImin3:ixImax3,iw)
   case ("ahotaka")
      ! First overwrite cells in domain with zero
      ixIpmin1=ixImin1;ixIpmin2=ixImax2+1;ixIpmin3=ixImin3;
      ixIpmax1=ixImax1;ixIpmax2=ixImax2+dixHotaka;ixIpmax3=ixImax3;
      do ix2=ixIpmin2,ixIpmax2
         w(ixIpmin1:ixIpmax1,ix2,ixIpmin3:ixIpmax3,iw) = zero
      end do
      ! Now apply asymm boundary
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         =-w(ixImin1:ixImax1,ixImax2+dixB:ixImax2+1:-1,ixImin3:ixImax3,iw)
   case ("special")
      ! skip it here, do AFTER all normal type boundaries are set
   case ("aperiodic")
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         = - w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw)
   case ("periodic")
      !            call mpistop("periodic bc info should come from neighbors")
   case default
      write (unitterm,*) "Undefined boundarytype ",typeB(iw,iB), &
           "for variable iw=",iw," and side iB=",iB
   end select


end do


      
      ! Apply on primitive variables?
      if (primitiveB(iside,idims)) then
         w(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) &
            = wsave(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) !restore from backup of the state in domain since p2c/c2p might lead to slight changes in w
         call conserve(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
            ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,w,x,patchfalse)
      end if
      
   end if 
case (3)
if (iside==2) then

   
   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! maximal boundary
   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      iB=ismax3
      ixImin1=ixBmin1;ixImin2=ixBmin2;ixImin3=ixBmax3+1-dixB;
      ixImax1=ixBmax1;ixImax2=ixBmax2;ixImax3=ixBmax3;

      !!!!!!! Make tangential ranges greater for staggered components
      !!!!!!! Add fill the normal component

      ! Apply on primitive variables?
      if (primitiveB(iside,idims)) then
         ixDmin1=ixBmin1;ixDmin2=ixBmin2;ixDmin3=ixBmin3+dixB;
         ixDmax1=ixBmax1;ixDmax2=ixBmax2;ixDmax3=ixBmax3-dixB;
         wsave(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) &
            = w(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) !make a backup of the state in domain since p2c/c2p might lead to slight changes in w
         call primitive(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
            ixDmin1,ixDmin2,ixDmin3,ixDmax1,ixDmax2,ixDmax3,w,x)
      end if
      
      ! cont/symm/asymm types
      do iw=1,nwflux+nwaux

         select case (typeB(iw,iB))
         case ("symm")
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               = w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3-1:ixImin3-dixB:-1,&
               iw)
         case ("asymm")
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               =-w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3-1:ixImin3-dixB:-1,&
               iw)
         case ("cont")
            do ix3=ixImin3,ixImax3
               w(ixImin1:ixImax1,ixImin2:ixImax2,ix3,iw) = w(ixImin1:ixImax1,&
                  ixImin2:ixImax2,ixImin3-1,iw)
            end do
         case("noinflow")
            if (iw==1+3)then
              do ix3=ixImin3,ixImax3
                  w(ixImin1:ixImax1,ixImin2:ixImax2,ix3,iw) &
                     = max(w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3-1,iw),&
                     zero)
              end do
            else
              do ix3=ixImin3,ixImax3
                  w(ixImin1:ixImax1,ixImin2:ixImax2,ix3,iw) &
                     = w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3-1,iw)
              end do
            end if
         case("limitinflow")
            if (iw==1+3)then
              do ix3=ixImin3,ixImax3
                  w(ixImin1:ixImax1,ixImin2:ixImax2,ix3,iw) &
                     = max(w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3-1,iw), &
                                           w(ixImin1:ixImax1,ixImin2:ixImax2,&
                                              ixImin3-1,iw)*ratebdflux)
              end do
            else
              do ix3=ixImin3,ixImax3
                  w(ixImin1:ixImax1,ixImin2:ixImax2,ix3,iw) &
                     = w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3-1,iw)
              end do
           end if
        case ("polefix")
           ! First overwrite cells in domain
           ixIpmin1=ixImin1;ixIpmin2=ixImin2;ixIpmin3=ixImin3-dixPolefix;
           ixIpmax1=ixImax1;ixIpmax2=ixImax2;ixIpmax3=ixImin3-1;
           do ix3=ixIpmin3,ixIpmax3
              w(ixIpmin1:ixIpmax1,ixIpmin2:ixIpmax2,ix3,iw) &
                 = w(ixIpmin1:ixIpmax1,ixIpmin2:ixIpmax2,ixIpmin3-1,iw)
           end do
           ! Now apply symm boundary
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               = w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3-1:ixImin3-dixB:-1,&
               iw)
        case ("apolefix")
           ! First interpolate cells in domain down to zero at pole
           ixIpmin1=ixImin1;ixIpmin2=ixImin2;ixIpmin3=ixImin3-dixPolefix;
           ixIpmax1=ixImax1;ixIpmax2=ixImax2;ixIpmax3=ixImax3-1;
           do ix3=ixIpmin3,ixIpmax3
              w(ixIpmin1:ixIpmax1,ixIpmin2:ixIpmax2,ix3,iw) &
                 = w(ixIpmin1:ixIpmax1,ixIpmin2:ixIpmax2,ixIpmin3-1,iw) &
                   * (xprobmax3-x(ixIpmin1:ixIpmax1,ixIpmin2:ixIpmax2,ix3,&
                      3)) &
                   / (xprobmax3-x(ixIpmin1:ixIpmax1,ixIpmin2:ixIpmax2,&
                      ixIpmin3-1,3))
           end do
        case ("hotaka")
           ! First overwrite cells in domain with zero
           ixIpmin1=ixImin1;ixIpmin2=ixImin2;ixIpmin3=ixImin3-dixHotaka;
           ixIpmax1=ixImax1;ixIpmax2=ixImax2;ixIpmax3=ixImin3-1;
           do ix3=ixIpmin3,ixIpmax3
              w(ixIpmin1:ixIpmax1,ixIpmin2:ixIpmax2,ix3,iw) = zero
           end do
           ! Now apply symm boundary
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               = w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3-1:ixImin3-dixB:-1,&
               iw)
        case ("ahotaka")
            ! First overwrite cells in domain with zero
           ixIpmin1=ixImin1;ixIpmin2=ixImin2;ixIpmin3=ixImin3-dixHotaka;
           ixIpmax1=ixImax1;ixIpmax2=ixImax2;ixIpmax3=ixImax3-1;
           do ix3=ixIpmin3,ixIpmax3
              w(ixIpmin1:ixIpmax1,ixIpmin2:ixIpmax2,ix3,iw) = zero
           end do
           ! Now apply asymm boundary
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               =-w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3-1:ixImin3-dixB:-1,&
               iw)
         case ("special")
            ! skip it here, do AFTER all normal type boundaries are set
         case ("aperiodic")
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
            w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
               = - w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw)
         case ("periodic")
!            call mpistop("periodic bc info should come from neighbors")
         case default
            write (unitterm,*) "Undefined boundarytype ",typeB(iw,iB), &
               "for variable iw=",iw," and side iB=",iB
         end select


      end do


      ! Apply on primitive variables?
      if (primitiveB(iside,idims)) then
         w(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) &
            = wsave(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) !restore from backup of the state in domain since p2c/c2p might lead to slight changes in w
         call conserve(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
            ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,w,x,patchfalse)
      end if
      
   else
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! minimal boundary
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      iB=ismin3
      ixImin1=ixBmin1;ixImin2=ixBmin2;ixImin3=ixBmin3;
      ixImax1=ixBmax1;ixImax2=ixBmax2;ixImax3=ixBmin3-1+dixB;

      ! Apply on primitive variables?
      if (primitiveB(iside,idims)) then         
         ixDmin1=ixBmin1;ixDmin2=ixBmin2;ixDmin3=ixBmin3+dixB;
         ixDmax1=ixBmax1;ixDmax2=ixBmax2;ixDmax3=ixBmax3-dixB;
         wsave(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) &
            = w(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) !make a backup of the state in domain since p2c/c2p might lead to slight changes in w
         call primitive(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
            ixDmin1,ixDmin2,ixDmin3,ixDmax1,ixDmax2,ixDmax3,w,x)
      end if
 
      ! cont/symm/asymm types
      do iw=1,nwflux+nwaux

   select case (typeB(iw,iB))
   case ("symm")
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         = w(ixImin1:ixImax1,ixImin2:ixImax2,ixImax3+dixB:ixImax3+1:-1,iw)
   case ("asymm")
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         =-w(ixImin1:ixImax1,ixImin2:ixImax2,ixImax3+dixB:ixImax3+1:-1,iw)
   case ("cont")
      do ix3=ixImin3,ixImax3
         w(ixImin1:ixImax1,ixImin2:ixImax2,ix3,iw) = w(ixImin1:ixImax1,&
            ixImin2:ixImax2,ixImax3+1,iw)
      end do
   case("noinflow")
      if (iw==1+3)then
         do ix3=ixImin3,ixImax3
            w(ixImin1:ixImax1,ixImin2:ixImax2,ix3,iw) = min(w(ixImin1:ixImax1,&
               ixImin2:ixImax2,ixImax3+1,iw),zero)
         end do
      else
         do ix3=ixImin3,ixImax3
            w(ixImin1:ixImax1,ixImin2:ixImax2,ix3,iw) = w(ixImin1:ixImax1,&
               ixImin2:ixImax2,ixImax3+1,iw)
         end do
      end if
   case("limitinflow")
      if (iw==1+3)then
         do ix3=ixImin3,ixImax3
            w(ixImin1:ixImax1,ixImin2:ixImax2,ix3,iw) = min(w(ixImin1:ixImax1,&
               ixImin2:ixImax2,ixImax3+1,iw), &
                 w(ixImin1:ixImax1,ixImin2:ixImax2,ixImax3+1,iw)*ratebdflux)
         end do
      else
         do ix3=ixImin3,ixImax3
            w(ixImin1:ixImax1,ixImin2:ixImax2,ix3,iw) = w(ixImin1:ixImax1,&
               ixImin2:ixImax2,ixImax3+1,iw)
         end do
      end if
   case ("polefix")
      ! First overwrite cells in domain
      ixIpmin1=ixImin1;ixIpmin2=ixImin2;ixIpmin3=ixImax3+1;
      ixIpmax1=ixImax1;ixIpmax2=ixImax2;ixIpmax3=ixImax3+dixPolefix;
      do ix3=ixIpmin3,ixIpmax3
         w(ixIpmin1:ixIpmax1,ixIpmin2:ixIpmax2,ix3,iw) = w(ixIpmin1:ixIpmax1,&
            ixIpmin2:ixIpmax2,ixIpmax3+1,iw)
      end do
      ! Now apply symm boundary
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         = w(ixImin1:ixImax1,ixImin2:ixImax2,ixImax3+dixB:ixImax3+1:-1,iw)
   case ("apolefix")
      ! First interpolate cells in domain down to zero at pole
      ixIpmin1=ixImin1;ixIpmin2=ixImin2;ixIpmin3=ixImax3+1;
      ixIpmax1=ixImax1;ixIpmax2=ixImax2;ixIpmax3=ixImax3+dixPolefix;
      do ix3=ixIpmin3,ixIpmax3
         w(ixIpmin1:ixIpmax1,ixIpmin2:ixIpmax2,ix3,iw) = w(ixIpmin1:ixIpmax1,&
            ixIpmin2:ixIpmax2,ixIpmax3+1,iw) &
              * (x(ixIpmin1:ixIpmax1,ixIpmin2:ixIpmax2,ix3,3))&
                 /x(ixIpmin1:ixIpmax1,ixIpmin2:ixIpmax2,ixIpmax3+1,3)
      end do
      ! Now apply asymm boundary
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         =-w(ixImin1:ixImax1,ixImin2:ixImax2,ixImax3+dixB:ixImax3+1:-1,iw)
   case ("hotaka")
      ! First overwrite cells in domain with zero
      ixIpmin1=ixImin1;ixIpmin2=ixImin2;ixIpmin3=ixImax3+1;
      ixIpmax1=ixImax1;ixIpmax2=ixImax2;ixIpmax3=ixImax3+dixHotaka;
      do ix3=ixIpmin3,ixIpmax3
         w(ixIpmin1:ixIpmax1,ixIpmin2:ixIpmax2,ix3,iw) = zero
      end do
      ! Now apply symm boundary
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         = w(ixImin1:ixImax1,ixImin2:ixImax2,ixImax3+dixB:ixImax3+1:-1,iw)
   case ("ahotaka")
      ! First overwrite cells in domain with zero
      ixIpmin1=ixImin1;ixIpmin2=ixImin2;ixIpmin3=ixImax3+1;
      ixIpmax1=ixImax1;ixIpmax2=ixImax2;ixIpmax3=ixImax3+dixHotaka;
      do ix3=ixIpmin3,ixIpmax3
         w(ixIpmin1:ixIpmax1,ixIpmin2:ixIpmax2,ix3,iw) = zero
      end do
      ! Now apply asymm boundary
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         =-w(ixImin1:ixImax1,ixImin2:ixImax2,ixImax3+dixB:ixImax3+1:-1,iw)
   case ("special")
      ! skip it here, do AFTER all normal type boundaries are set
   case ("aperiodic")
 !this just multiplies the variables with (-), they have been set from neighbors just like periodic.
      w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw) &
         = - w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iw)
   case ("periodic")
      !            call mpistop("periodic bc info should come from neighbors")
   case default
      write (unitterm,*) "Undefined boundarytype ",typeB(iw,iB), &
           "for variable iw=",iw," and side iB=",iB
   end select


end do


      
      ! Apply on primitive variables?
      if (primitiveB(iside,idims)) then
         w(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) &
            = wsave(ixDmin1:ixDmax1,ixDmin2:ixDmax2,ixDmin3:ixDmax3,1:nw) !restore from backup of the state in domain since p2c/c2p might lead to slight changes in w
         call conserve(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
            ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,w,x,patchfalse)
      end if
      
   end if 
end select


! do special case AFTER all normal cases are set
if (any(typeB(1:nwflux+nwaux,iB)=="special")) then
   call specialbound_usr(time,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
      ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,iB,s)
end if

   end associate
end subroutine bc_phys
!=============================================================================
subroutine getintbc(time,psuse)

use mod_amrvacdef

double precision, intent(in)              :: time
type(state), dimension(ngridshi)          :: psuse

! .. local ..
integer :: iigrid, igrid, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,&
   level
!----------------------------------------------------------------------------
ixOmin1=ixGlo1+dixB;ixOmin2=ixGlo2+dixB;ixOmin3=ixGlo3+dixB
ixOmax1=ixGhi1-dixB;ixOmax2=ixGhi2-dixB;ixOmax3=ixGhi3-dixB;

!$OMP PARALLEL DO SCHEDULE(RUNTIME) PRIVATE(igrid,level)
do iigrid=1,igridstail_active; igrid=igrids_active(iigrid);
   call set_tmpGlobals(igrid)
   level = node(plevel_,igrid)
   call bc_int(level,time,ixGlo1,ixGlo2,ixGlo3,ixGhi1,ixGhi2,ixGhi3,ixOmin1,&
      ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,psuse(igrid)%w%w,px(igrid)%x)
end do
!$OMP END PARALLEL DO

      
end subroutine getintbc
!=============================================================================
subroutine is_neighbor_coarse(s,is_coarse)

use mod_amrvacdef

type(state), intent(in)          :: s
logical, intent(out)             :: is_coarse(ndim,2)
! .. local ..
integer                          :: i1,i2,i3, i
!-----------------------------------------------------------------------------

is_coarse=.false.
! If we are the coarse version of a buffer
! dont consider neighbor as coarser
if (s%is_coarse .eqv. .true.) return 

   do i3=-1,1
   do i2=-1,1
   do i1=-1,1
      if ((i1.eq.0.and.i2.eq.0.and.i3.eq.0).or.(abs(i1)+abs(i2)+abs&
         (i3)).ne.1) cycle
      if (neighbor_type(i1,i2,i3,s%igrid).eq.2) then
         if (i1.eq.-1) is_coarse(1,1) = .true.
         if (i2.eq.-1) is_coarse(2,1) = .true.
         if (i3.eq.-1) is_coarse(3,1) = .true.
         if (i1.eq.+1) is_coarse(1,2) = .true.
         if (i2.eq.+1) is_coarse(2,2) = .true.
         if (i3.eq.+1) is_coarse(3,2) = .true.
      end if
   end do
   end do
   end do

end subroutine is_neighbor_coarse
!=============================================================================
