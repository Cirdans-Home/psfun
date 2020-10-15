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
module psfun_d_krylov_mod
    ! The psfun_d_krylov_mod contains the generic call to a Krylov subspace
    ! method for the computation of :math:`y = f(A) x`, for :math:`A` large and
    ! sparse.

  use psb_base_mod
  use psfun_d_serial_mod
  implicit none

  type, public :: psfun_d_krylov
    character(len=20)   :: kname   = 'ARNOLDI' ! Name of the Krylov method
  contains
    ! Set the options
    procedure, pass(fun) :: setstring  => psfun_d_setstring
    generic, public :: set => setstring
  end type

  private :: psfun_d_setstring

contains

  subroutine psfun_d_setstring(fun,what,val,info)
      ! Set function for setting options defined by a string
      use psb_base_mod
      implicit none

      class(psfun_d_krylov), intent(inout) :: fun   ! Function object
      character(len=*), intent(in)         :: what  ! String of option to set
      character(len=*), intent(in)         :: val   ! Value of the string
      integer(psb_ipk_), intent(out)       :: info  ! Output flag

      info = psb_success_
      select case (psb_toupper(what))
      case ("KNAME")
          fun%kname = val
      case default
          info = psb_err_invalid_args_combination_
      end select

  end subroutine psfun_d_setstring

end module psfun_d_krylov_mod
