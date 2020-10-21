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
module psfun_z_quadrature_mod
! This module computes the matrix-function vector product by means of the
! approximation of :math:`f(A)\mathbf{x}` based on quadrature formula, i.e.,
! having computed the poles and the scalings of the formula solves N linear
! systems to approximate the product.

  use psfun_base_quadrature_mod
  use psb_base_mod
  use psfun_utils_mod
  implicit none

  type, extends(psfun_quadrature) :: psfun_z_quadrature
    complex(psb_dpk_), allocatable, dimension(:) :: xi     ! Poles of the formula
    complex(psb_dpk_), allocatable, dimension(:) :: c      ! Scaling of the formula
    real(psb_dpk_)                               :: eta    ! Global Scaling
    real(psb_dpk_)                               :: sign   ! Sign for A
    type(psb_dspmat_type), pointer               :: a      ! Matrix on which we work
    type(psb_dprec_type), pointer                :: prec   ! Preconditioner for the solution of the associate linear systems
  contains
    procedure, pass(quad) :: computepoles       => psfun_z_computepoles
    procedure, pass(quad) :: setmatrix          => psfun_z_setmatrix
    procedure, pass(quad) :: setpreconditioner  => psfun_z_setpreconditioner
    generic, public :: set => setmatrix, setpreconditioner
  end type psfun_z_quadrature

  ! ************************************************************************** !
  ! Abstract interfaces
  ! These are abstract interfaces that can be used to implement any possible
  ! quadrature rule on the user side. Then the computepoles routine takes
  ! it as input together with the number of wanted poles to populate the
  ! the variable of type psfun_quadrature.
  abstract interface
    subroutine zquadrule(zfun,xi,c,eta,sign,N,info,cparams,rparams)
      ! To integrate a function that take as inputs complex number and gives as output
      ! complex numbers
      use psb_base_mod
      implicit none
      complex(psb_dpk_), allocatable, dimension(:), intent(out) :: xi     ! Poles of the formula
      complex(psb_dpk_), allocatable, dimension(:), intent(out) :: c      ! Scaling of the formula
      real(psb_dpk_), intent(out)               :: eta! Global Scaling
      real(psb_dpk_), intent(out)               :: sign   ! Sign for A
      procedure (zquadfun), pointer, intent(in) :: zfun    ! Function to integrate
      integer(psb_ipk_), intent(in)             :: N ! Number of Poles
      integer(psb_ipk_), intent(out)            :: info ! Flag on the results
      complex(psb_dpk_), dimension(:), optional, intent(in) :: cparams ! Optional complex parameters
      real(psb_dpk_), dimension(:), optional, intent(in)    :: rparams ! Optional real parameters
    end subroutine
  end interface
  ! ************************************************************************** !
  ! These re abstract interfaces for the function to which the quadrature rule
  ! should be applied
  abstract interface
    function zquadfun(z) result(res)
      use psb_base_mod
      implicit none
      complex(psb_dpk_), intent(in) :: z
      complex(psb_dpk_) :: res
    end function
  end interface
  ! ************************************************************************** !

  ! Module procedures implementing different quadrature rules, all the rules
  ! are contained in the relative submodule
  interface hhtmethod1
    module subroutine hhtmethod1(zfun,xi,c,eta,sign,N,info,rparams)
      procedure (zquadfun), pointer, intent(in)                 :: zfun   ! Function to integrate
      complex(psb_dpk_), allocatable, dimension(:), intent(out) :: xi     ! Poles of the formula
      complex(psb_dpk_), allocatable, dimension(:), intent(out) :: c      ! Scaling of the formula
      real(psb_dpk_), intent(out)     :: eta! Global Scaling
      real(psb_dpk_), intent(out)     :: sign   ! Sign for A
      integer(psb_ipk_), intent(in)   :: N ! Number of Poles
      integer(psb_ipk_), intent(out)  :: info ! Flag on the results
      real(psb_dpk_), dimension(2), intent(in) :: rparams ! Optional real parameters
    end subroutine
  end interface

contains

  subroutine psfun_z_computepoles(quad,quadformula,zfun,N,info,cparams,rparams)
    use psb_base_mod
    implicit none

    class(psfun_z_quadrature), intent(inout)    :: quad      ! Quadrature type
    procedure (zquadrule), pointer, intent(in):: quadformula ! Quadrature formula
    procedure (zquadfun), pointer, intent(in) :: zfun    ! Function to integrate
    integer(psb_ipk_), intent(in)             :: N           ! Number of poles
    integer(psb_ipk_), intent(out)            :: info        ! Flag on the results
    complex(psb_dpk_), dimension(:), optional, intent(in) :: cparams ! Optional complex parameters
    real(psb_dpk_), dimension(:), optional, intent(in) :: rparams ! Optional real parameters

    info = 0

    if( present(cparams) ) then
      call quadformula(zfun,quad%xi,quad%c,quad%eta,quad%sign,N,info,cparams=cparams)
    else if( present(rparams) ) then
      call quadformula(zfun,quad%xi,quad%c,quad%eta,quad%sign,N,info,rparams=rparams)
    else if( (present(cparams)).and.(present(rparams)) ) then
      call quadformula(zfun,quad%xi,quad%c,quad%eta,quad%sign,N,info,cparams,rparams)
    else
      call quadformula(zfun,quad%xi,quad%c,quad%eta,quad%sign,N,info)
    end if

    return

  end subroutine psfun_z_computepoles

  subroutine psfun_z_setmatrix(quad,a)
    use psb_base_mod
    implicit none

    class(psfun_z_quadrature), intent(inout)    :: quad   ! Quadrature type
    type(psb_dspmat_type), target               :: a      ! Matrix on which we work

    quad%a => a

  end subroutine psfun_z_setmatrix

  subroutine psfun_z_setpreconditioner(quad,prec)
    use psb_base_mod
    implicit none

    class(psfun_z_quadrature), intent(inout)    :: quad   ! Quadrature type
    type(psb_dprec_type), target                :: prec   ! Preconditioner for the solution of the associate linear systems

    quad%prec => prec

  end subroutine psfun_z_setpreconditioner

end module psfun_z_quadrature_mod
