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

module mod_amr_fct
  use mod_indices, only: ngridshi
  implicit none
  
  type facealloc
    double precision, dimension(:,:,:), pointer :: face
  end type facealloc

  type fake_neighbors
    integer :: igrid
    integer :: ipe
  end type fake_neighbors

  type(facealloc), dimension(2,3,ngridshi), save :: pface

  type(fake_neighbors), dimension(2,2,2,3,ngridshi), save :: fine_neighbors

  integer, dimension(2,-1:11,-1:12,-1:13,ngridshi), save :: old_neighbor

  integer, save        :: ibuf_recv, ibuf_send, ibuf_send_next

contains


end module mod_amr_fct

