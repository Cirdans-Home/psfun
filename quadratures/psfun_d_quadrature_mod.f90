! BSD 3-Clause License
!
! Copyright (c) 2020, Fabio Durastante
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this
!    list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its
!    contributors may be used to endorse or promote products derived from
!    this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
module psfun_d_quadrature_mod
! This module computes the matrix-function vector product by means of the
! approximation of :math:`f(A)\mathbf{x}` based on quadrature formula, i.e.,
! having computed the poles and the scalings of the formula solves N linear
! systems to approximate the product.

  use psfun_base_quadrature_mod
  use psb_base_mod
  use psfun_utils_mod
  implicit none

  type, extends(psfun_quadrature) :: psfun_d_quadrature
    real(psb_dpk_), allocatable, dimension(:) :: xi     ! Poles of the formula
    real(psb_dpk_), allocatable, dimension(:) :: c      ! Scaling of the formula
    real(psb_dpk_)                            :: eta    ! Global Scaling
    real(psb_dpk_)                            :: sign   ! Sign for A
    type(psb_dspmat_type), pointer :: a    ! Matrix on which we work
    type(psb_dprec_type)           :: prec ! Preconditioner for the solution of the associated linear systems
  contains

  end type psfun_d_quadrature

  abstract interface
    subroutine dquadrule(dfun,xi,c,eta,sign,N,info,cparams,rparams)
      ! To integrate a function that take as inputs real numbers and gives as
      ! output complex numbers
      use psb_base_mod
      implicit none
      complex(psb_dpk_), allocatable, dimension(:), intent(out) :: xi     ! Poles of the formula
      complex(psb_dpk_), allocatable, dimension(:), intent(out) :: c      ! Scaling of the formula
      complex(psb_dpk_), intent(out) :: eta! Global Scaling
      complex(psb_dpk_), intent(out) :: sign   ! Sign for A
      procedure (dquadfun), pointer, intent(in) :: dfun    ! Function to integrate
      integer(psb_ipk_), intent(in)  :: N ! Number of Poles
      integer(psb_ipk_), intent(out) :: info ! Flag on the results
      complex(psb_dpk_), dimension(:), optional, intent(in) :: cparams ! Optional complex parameters
      real(psb_dpk_), dimension(:), optional, intent(in) :: rparams ! Optional real parameters
    end subroutine
  end interface

  abstract interface
    function dquadfun(z) result(res)
      use psb_base_mod
      implicit none
      real(psb_dpk_), intent(in) :: z
      complex(psb_dpk_) :: res
    end function
  end interface

contains

end module psfun_d_quadrature_mod
