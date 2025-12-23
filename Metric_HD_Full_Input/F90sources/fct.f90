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


subroutine fct_average(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
    ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, fC)

! Performs the average for the flux constrained transport from Toth 2000. 
! His equation (25)

use mod_amrvacdef

integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision, intent(inout)    :: fC(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nwflux,1:ndim)

double precision                   :: fB(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndir,1:ndim)
integer                            :: iwdim, iwdir, idim, idir
integer                            :: ixJpmin1,ixJpmin2,ixJpmin3,ixJpmax1,&
   ixJpmax2,ixJpmax3, ixKpmin1,ixKpmin2,ixKpmin3,ixKpmax1,ixKpmax2,ixKpmax3,&
    ixKmmin1,ixKmmin2,ixKmmin3,ixKmmax1,ixKmmax2,ixKmmax3, ixCmin1,ixCmin2,&
   ixCmin3,ixCmax1,ixCmax2,ixCmax3, ixJpKmmin1,ixJpKmmin2,ixJpKmmin3,&
   ixJpKmmax1,ixJpKmmax2,ixJpKmmax3
!-----------------------------------------------------------------------------


do idim = 1, ndim ! idim is the direction of the flux interface
   do idir = 1, ndim ! idir is the component of the field to consider
      if (idim.eq.idir) cycle
      iwdir = b0_+idir; iwdim = b0_+idim

! Assemble indices:
      ixCmax1=ixOmax1+1;ixCmax2=ixOmax2+1;ixCmax3=ixOmax3+1;
      ixCmin1=ixOmin1-1-kr(idim,1);ixCmin2=ixOmin2-1-kr(idim,2)
      ixCmin3=ixOmin3-1-kr(idim,3); !Extend range in flux direction
 
      ixJpmin1=ixCmin1+kr(idim,1);ixJpmin2=ixCmin2+kr(idim,2)
      ixJpmin3=ixCmin3+kr(idim,3);ixJpmax1=ixCmax1+kr(idim,1)
      ixJpmax2=ixCmax2+kr(idim,2);ixJpmax3=ixCmax3+kr(idim,3);
      ixKpmin1=ixCmin1+kr(idir,1);ixKpmin2=ixCmin2+kr(idir,2)
      ixKpmin3=ixCmin3+kr(idir,3);ixKpmax1=ixCmax1+kr(idir,1)
      ixKpmax2=ixCmax2+kr(idir,2);ixKpmax3=ixCmax3+kr(idir,3);
      ixKmmin1=ixCmin1-kr(idir,1);ixKmmin2=ixCmin2-kr(idir,2)
      ixKmmin3=ixCmin3-kr(idir,3);ixKmmax1=ixCmax1-kr(idir,1)
      ixKmmax2=ixCmax2-kr(idir,2);ixKmmax3=ixCmax3-kr(idir,3);
      ixJpKmmin1=ixJpmin1-kr(idir,1);ixJpKmmin2=ixJpmin2-kr(idir,2)
      ixJpKmmin3=ixJpmin3-kr(idir,3);ixJpKmmax1=ixJpmax1-kr(idir,1)
      ixJpKmmax2=ixJpmax2-kr(idir,2);ixJpKmmax3=ixJpmax3-kr(idir,3);


! Perform flux average:
         fB(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir,idim) = &
              0.125d0 * ( 2.0d0* fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
                 ixCmin3:ixCmax3,iwdir,idim) &
              + fC(ixKpmin1:ixKpmax1,ixKpmin2:ixKpmax2,ixKpmin3:ixKpmax3,&
                 iwdir,idim) &
              + fC(ixKmmin1:ixKmmax1,ixKmmin2:ixKmmax2,ixKmmin3:ixKmmax3,&
                 iwdir,idim) - dxlevel(idir) &
              * (fC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,iwdim,&
                 idir)  &
              + fC(ixJpmin1:ixJpmax1,ixJpmin2:ixJpmax2,ixJpmin3:ixJpmax3,&
                 iwdim,idir)  &
              + fC(ixKmmin1:ixKmmax1,ixKmmin2:ixKmmax2,ixKmmin3:ixKmmax3,&
                 iwdim,idir)  &
              + fC(ixJpKmmin1:ixJpKmmax1,ixJpKmmin2:ixJpKmmax2,&
                 ixJpKmmin3:ixJpKmmax3,iwdim,idir)  &
              )/dxlevel(idim))

   end do ! idir
end do ! idim

! Overwrite the new flux entries:
do idim = 1, ndim ! idim is the direction of the flux interface
   do idir = 1, ndim ! idir is the component of the field to consider
      iwdir = b0_+idir
      if (idim.eq.idir) then
         ! To conserve divb to machine precision, this needs to be zero.
         ! However, since restriction and prolongation can introduce divb
 !when using AMR, additional measures for divb-control must be taken. 
         ! These rely on the normal field flux component (e.g. as in GLM).
         ! Thus when more than one level is used, we dont set this zero 
         ! to allow additional divb-control to work.
         if (mxnest .eq. 1) fC(ixImin1:ixImax1,ixImin2:ixImax2,&
            ixImin3:ixImax3,iwdir,idim) = zero
         cycle
      end if

      fC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,iwdir,idim) &
         = fB(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,idir,idim)

   end do ! idir
end do ! idim

end subroutine fct_average

!=============================================================================

!=============================================================================

subroutine b_from_vectorpotential(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, w, x)

use mod_amrvacdef

integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3, ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision, intent(inout)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:nw)
double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndim)

integer                            :: ixCmin1,ixCmin2,ixCmin3,ixCmax1,ixCmax2,&
   ixCmax3, ixCpmin1,ixCpmin2,ixCpmin3,ixCpmax1,ixCpmax2,ixCpmax3, ixCmmin1,&
   ixCmmin2,ixCmmin3,ixCmmax1,ixCmmax2,ixCmmax3, ixOmmin1,ixOmmin2,ixOmmin3,&
   ixOmmax1,ixOmmax2,ixOmmax3, hxOmin1,hxOmin2,hxOmin3,hxOmax1,hxOmax2,&
   hxOmax3, idim, idir, j, k
double precision                   :: xC(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:ndim), A(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,&
   1:ndir), dxC(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,1:ndim)
double precision                   :: B(ixImin1:ixImax1,ixImin2:ixImax2,&
   ixImin3:ixImax3,1:3)
!-----------------------------------------------------------------------------

A(:,:,:,:)=zero

ixCmax1=ixOmax1;ixCmax2=ixOmax2;ixCmax3=ixOmax3;
ixCmin1=ixOmin1-1;ixCmin2=ixOmin2-1;ixCmin3=ixOmin3-1; ! Extend range by one

do idir=7-2*ndim,ndir
 do idim=1,ndim
  ! Get corner coordinates
  if (idim.ne.idir) then
   ixCpmin1=ixCmin1+kr(idim,1);ixCpmin2=ixCmin2+kr(idim,2)
   ixCpmin3=ixCmin3+kr(idim,3);ixCpmax1=ixCmax1+kr(idim,1)
   ixCpmax2=ixCmax2+kr(idim,2);ixCpmax3=ixCmax3+kr(idim,3);
   ixCmmin1=ixCmin1-kr(idim,1);ixCmmin2=ixCmin2-kr(idim,2)
   ixCmmin3=ixCmin3-kr(idim,3);ixCmmax1=ixCmax1-kr(idim,1)
   ixCmmax2=ixCmax2-kr(idim,2);ixCmmax3=ixCmax3-kr(idim,3);
   xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim) &
      = half * (x(ixCpmin1:ixCpmax1,ixCpmin2:ixCpmax2,ixCpmin3:ixCpmax3,&
      idim) + x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim))
   dxC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim) &
      = xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
      idim) - xC(ixCmmin1:ixCmmax1,ixCmmin2:ixCmmax2,ixCmmin3:ixCmmax3,idim)
  else
   xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim)&
      =x(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idim)
  end if
 end do
 ! Initialise vector potencial at the corners

 call initvecpot_usr(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3, ixCmin1,&
    ixCmin2,ixCmin3,ixCmax1,ixCmax2,ixCmax3, xC, A(ixImin1:ixImax1,&
    ixImin2:ixImax2,ixImin3:ixImax3,idir), idir)

end do

! ---
! take the curl of the vectorpotential: 

B(:,:,:,:) = zero
do idir = 1, 3
   do j = 1, ndim
      if (idir .eq. j) cycle
      ixCmin1=ixOmin1-kr(idir,1);ixCmin2=ixOmin2-kr(idir,2)
      ixCmin3=ixOmin3-kr(idir,3);
      ixCmmin1=ixCmin1-kr(j,1);ixCmmin2=ixCmin2-kr(j,2)
      ixCmmin3=ixCmin3-kr(j,3);ixCmmax1=ixCmax1-kr(j,1)
      ixCmmax2=ixCmax2-kr(j,2);ixCmmax3=ixCmax3-kr(j,3);
      do k = 1,ndir
         ! Field on the faces
         select case(idir)
         case(1)
         
         
         where (mygeo%surfaceC1(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
            ixCmin3:ixCmax3) .ne. zero)
            B(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir) &
               = B(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir) &
                 + lvc(idir,j,k) * dxlevel(k) &
                 * (A(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                    k)-A(ixCmmin1:ixCmmax1,ixCmmin2:ixCmmax2,&
                    ixCmmin3:ixCmmax3,k)) 
         elsewhere
            B(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir) = zero
         end where
        
         
case(2)
         
         
         where (mygeo%surfaceC2(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
            ixCmin3:ixCmax3) .ne. zero)
            B(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir) &
               = B(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir) &
                 + lvc(idir,j,k) * dxlevel(k) &
                 * (A(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                    k)-A(ixCmmin1:ixCmmax1,ixCmmin2:ixCmmax2,&
                    ixCmmin3:ixCmmax3,k)) 
         elsewhere
            B(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir) = zero
         end where
        
         
case(3)
         
         
         where (mygeo%surfaceC3(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
            ixCmin3:ixCmax3) .ne. zero)
            B(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir) &
               = B(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir) &
                 + lvc(idir,j,k) * dxlevel(k) &
                 * (A(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,&
                    k)-A(ixCmmin1:ixCmmax1,ixCmmin2:ixCmmax2,&
                    ixCmmin3:ixCmmax3,k)) 
         elsewhere
            B(ixCmin1:ixCmax1,ixCmin2:ixCmax2,ixCmin3:ixCmax3,idir) = zero
         end where
        
         
         end select
      end do
   end do
end do



! Average to the cell centers and fill solution array:

do idir = 1, 3
   ixOmmin1=ixOmin1-kr(idir,1);ixOmmin2=ixOmin2-kr(idir,2)
   ixOmmin3=ixOmin3-kr(idir,3);ixOmmax1=ixOmax1-kr(idir,1)
   ixOmmax2=ixOmax2-kr(idir,2);ixOmmax3=ixOmax3-kr(idir,3);
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,b0_+idir) &
      = half * dxC(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
      idir)*(B(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
      idir) + B(ixOmmin1:ixOmmax1,ixOmmin2:ixOmmax2,ixOmmin3:ixOmmax3,idir))&
      /mygeo%dvolume(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
end do

end subroutine b_from_vectorpotential
!=============================================================================


subroutine printarray(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,array)
! Routine for easily printing an array, just for debugging, erase later

use mod_amrvacdef

integer, intent(in)           :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
   ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
double precision, intent(in)  :: array(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
   ixOmin3:ixOmax3)
!double precision, intent(in)  :: array(:,:)
integer                       :: i,j,k
!----------------------------------------------------------------------------

print *, array


   do i=ixOmin1,ixOmax1




      do k=ixOmin3,ixOmax3
        write (*,'(I2,I2,I2,E17.9,$)') i,j,k,array(i,j,k)
      end do
      print *, ' '
     
   end do

   print *, '---'


write (*,'(A)') '------------------------------------------------------------'
end subroutine printarray
!=============================================================================
