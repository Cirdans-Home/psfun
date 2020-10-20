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
submodule (psfun_utils_mod) psfun_z_utils_mod
  use psb_base_mod
  implicit none

  complex, parameter :: iunit = cmplx(0.0_psb_dpk_,1.0_psb_dpk_)

contains

  recursive module subroutine z_ellipj(u,L,sn,cn,dn,flag)
    ! Returns the values of the Jacobi elliptic functions evaluated at double
    ! argument u and parameter :math:`M = \exp(-2 \pi L)`, :math:`0 < L < \infty`.
    ! For :math:`M = K^2`, and `K` the elliptic modulus.
    !
    complex(psb_dpk_), intent(in)    :: u
    real(psb_dpk_), intent(in)       :: L
    complex(psb_dpk_), intent(out)   :: sn,cn,dn
    logical, optional, intent(in)    :: flag

    real(psb_dpk_)      :: m,K(2),kappa,mu
    complex(psb_dpk_)   :: ucp, sinu, cosu, v, sn1, cn1, dn1, denom, snh, cnh, dnh
    logical             :: flag_, high
    real(psb_dpk_), parameter :: coeffs(7) = (/0.0,1.0,2.0,5.0,14.0,42.0,132.0/)

    ucp = u

    if( present(flag) ) then
      flag_ = flag
    else
      flag_ = .true.
    end if

    if(flag_) then
      K = ellipkkp(L)
      high = (aimag(ucp) > K(2)/2.0_psb_dpk_ )
      if (high) then
        ucp = iunit*K(2) - ucp
      end if
      m = exp(-2.0_psb_dpk_*DPI*L)
    else
      ! In case of recursive call we have already transformed L into m
      high = .false.
      m = L
    endif

    if (m < 4.0_psb_dpk_*EPSILON(m) ) then
      sinu = sin(u);
      cosu = cos(u);
      sn = sinu + m/4.0_psb_dpk_*( sinu*cosu-u)*cosu;
      cn = cosu + m/4.0_psb_dpk_*(-sinu*cosu+u)*sinu;
      dn = 1    + m/4.0_psb_dpk_*( cosu**2 -sinu**2-1.0_psb_dpk_);
    else
      if( m > 1e-3) then
        kappa = (1.0_psb_dpk_-sqrt(1.0_psb_dpk_-m))/(1.0_psb_dpk_+sqrt(1.0_psb_dpk_-m))
      else
        kappa = horner( coeffs , m/4.0_psb_dpk_ )
      end if
      mu = kappa**2
      v = u/(1.0_psb_dpk_ + kappa)
      call z_ellipj(v,mu,sn1,cn1,dn1,.false.)
      denom = (1.0_psb_dpk_ + kappa*sn1**2)
      sn = (1+kappa)*sn1/denom
      cn = cn1*dn1/denom
      dn = (1-kappa*sn1**2)/denom
    end if

    if(high) then
      snh = sn
      cnh = cn
      dnh = dn
      sn = -1.0_psb_dpk_/(sqrt(m)*snh)
      cn = iunit*dnh/(sqrt(m)*snh)
      dn = iunit*cnh/snh
    endif

    return

  end subroutine

end submodule
