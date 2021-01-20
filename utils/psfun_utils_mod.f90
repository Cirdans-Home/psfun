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
module psfun_utils_mod
  use psb_base_mod
  implicit none

  ! Fixed parameters
  real(psb_dpk_), private, parameter :: DPI = 4.0_psb_dpk_*ATAN(1.0_psb_dpk_) ! Double precision :math:`\pi` for internal usage

  interface ellipkkp
      module function ellipkkp(L) result(K)
        real(psb_dpk_), intent(in)  :: L
        real(psb_dpk_)              :: K(2)
      end function
  end interface ellipkkp

  interface ellipj
    module subroutine d_ellipj(u,L,sn,cn,dn)
      real(psb_dpk_), intent(in)    :: u
      real(psb_dpk_), intent(in)    :: L
      real(psb_dpk_), intent(out)   :: sn,cn,dn
    end subroutine
    recursive module subroutine z_ellipj(u,L,sn,cn,dn,flag)
      complex(psb_dpk_), intent(in)    :: u
      real(psb_dpk_), intent(in)       :: L
      complex(psb_dpk_), intent(out)   :: sn,cn,dn
      logical, optional, intent(in)    :: flag
    end subroutine
  end interface ellipj

  interface horner
    module function horner(coeffs, x) result (res)
      real(psb_dpk_), dimension (:), intent (in) :: coeffs
      real(psb_dpk_), intent (in) :: x
      real(psb_dpk_) :: res
    end function
  end interface horner

  interface linspace
    module subroutine linspace(from, to, array, range)
      real(psb_dpk_), intent(in) :: from, to
      real(psb_dpk_), intent(out) :: array(:)
      real(psb_dpk_), optional, intent(in) :: range
    end subroutine
  end interface

contains

end module psfun_utils_mod
