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
submodule (psfun_utils_mod) psfun_d_utils_mod
  !! Real variants of the utils functions
  use psb_base_mod
  implicit none

contains

  module function ellipkkp(L) result(K)
    !! Complete elliptic integral of the first kind, with complement.
    !! Returns the value of the complete elliptic integral of the first kind,
    !! evaluated at \(M=\exp(-2\pi L)\), \(0 < L < \infty\), and the
    !! complementarity parameter \(1-M\).
    use psb_base_mod
    implicit none

    real(psb_dpk_), intent(in)  :: L
    real(psb_dpk_)              :: K(2)

    ! Local variables
    real(psb_dpk_)    :: M,a0,a1,b0,b1,c1,w1,s0,MM
    integer(psb_ipk_) :: i1

    ! When M=exp(-2 π L) ≈ 0, use an O(M) approximation
    if(L > 10) then
      K(1) = DPI/2_psb_dpk_
      K(2) = DPI*L + log(4.0_psb_dpk_)
      return
    end if

    M = exp(-2.0_psb_dpk_*DPI*L)
    a0 = 1.0_psb_dpk_
    b0 = sqrt(1-M)
    s0 = M
    i1 = 0
    MM = 1.0_psb_dpk_
    do while (MM > EPSILON(MM) )
      a1 = (a0+b0)/2.0_psb_dpk_
      b1 = sqrt(a0*b0)
      c1 = (a0-b0)/2.0_psb_dpk_
      i1 = i1 + 1
      w1 = (2.0_psb_dpk_**i1)*(c1**2)
      MM = w1
      s0 = s0 + w1
      a0 = a1
      b0 = b1
    end do
    K(1) = DPI/(2.0_psb_dpk_*a1)

    a0 = 1.0_psb_dpk_
    b0 = sqrt(M)
    s0 = 1-M
    i1 = 0
    MM = 1.0_psb_dpk_
    do while (MM > EPSILON(MM) )
      a1 = (a0+b0)/2.0_psb_dpk_
      b1 = sqrt(a0*b0)
      c1 = (a0-b0)/2.0_psb_dpk_
      i1 = i1 + 1
      w1 = (2.0_psb_dpk_**i1)*(c1**2)
      MM = w1
      s0 = s0 + w1
      a0 = a1
      b0 = b1
    end do
    K(2) = DPI/(2.0_psb_dpk_*a1)


    return

  end function ellipkkp

  module subroutine d_ellipj(u,L,sn,cn,dn)
    !! Returns the values of the Jacobi elliptic functions evaluated at double
    !! argument u and parameter \(M = \exp(-2 \pi L)\), \(0 < L < \infty\).
    !! For \(M = K^2\), and \(K\) the elliptic modulus.
    !
    real(psb_dpk_), intent(in)    :: u
    real(psb_dpk_), intent(in)    :: L
    real(psb_dpk_), intent(out)   :: sn,cn,dn

    real(psb_dpk_) :: m
    m = exp(-2.0_psb_dpk_*DPI*L)

    call sncndn( u, m, sn, cn, dn )

    return

  end subroutine

  module function horner(coeffs, x) result (res)
    !! Apply Horner rule to evaluate a polynomial

    use psb_base_mod
    implicit none
    real(psb_dpk_), dimension (:), intent (in) :: coeffs
    real(psb_dpk_), intent (in) :: x
    real(psb_dpk_) :: res

    ! Local Variable
    integer :: i

    res = 0.0_psb_dpk_
    do i = size(coeffs), 1, -1
      res = res * x + coeffs (i)
    end do

    return

  end function horner

end submodule
