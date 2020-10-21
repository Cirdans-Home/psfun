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
submodule (psfun_z_quadrature_mod) psfun_z_quadrules_mod
  ! This submodule contains the impliementation of the legacy quadrature rules
  ! distributed with the library.

  ! Constants used internally
  real(psb_dpk_), parameter :: DPI = 4.0_psb_dpk_*ATAN(1.0_psb_dpk_) ! Double precision :math:`\pi` for internal usage
  complex, parameter :: iunit = cmplx(0.0_psb_dpk_,1.0_psb_dpk_) ! Immaginary unit

contains

  module subroutine hhtmethod1(zfun,xi,c,eta,sign,N,info,rparams)
    ! Method 1 of Hale, Nicholas; Higham, Nicholas J.; Trefethen, Lloyd N.
    ! Computing ${\bf A}^\alpha,\ \log({\bf A})$, and related matrix functions
    ! by contour integrals. SIAM J. Numer. Anal. 46 (2008), no. 5, 2505--2523.
    use psb_base_mod
    use psfun_utils_mod, only: ellipkkp, ellipj
    implicit none
    procedure (zquadfun), pointer, intent(in) :: zfun               ! Function to integrate
    complex(psb_dpk_), allocatable, dimension(:), intent(out) :: xi     ! Poles of the formula
    complex(psb_dpk_), allocatable, dimension(:), intent(out) :: c      ! Scaling of the formula
    real(psb_dpk_), intent(out)    :: eta! Global Scaling
    real(psb_dpk_), intent(out)    :: sign   ! Sign for A
    integer(psb_ipk_), intent(in)  :: N ! Number of Poles
    integer(psb_ipk_), intent(out) :: info ! Flag on the results
    real(psb_dpk_), dimension(2), intent(in) :: rparams ! Optional real parameters

    ! local variables
    real(psb_dpk_)    :: mineig,maxeig,cond,L,K,Kp,Kvec(2)
    complex(psb_dpk_) :: u,cn,dn
    complex(psb_dpk_) :: t
    integer(psb_ipk_) :: j

    info = 0

    mineig = min(rparams(1),rparams(2)) ! Estimate of the minimum eigenvalue
    maxeig = max(rparams(1),rparams(2)) ! Estimate of the maximum eigenvalue

    allocate(xi(N), stat=info)
    allocate(c(N), stat=info)

    cond = (sqrt(maxeig/mineig)-1)/(sqrt(maxeig/mineig)+1)
    L = -log(k)/DPI
    Kvec = ellipkkp(L)
    K = Kvec(1)
    Kp = Kvec(2)
    do j=1,N,1
      t = .5_psb_dpk_*iunit*Kp - K + (2.0_psb_dpk_*(j-1_psb_ipk_)+1_psb_dpk_)*K/N
      call ellipj(t,L,u,cn,dn)
      xi(j) = sqrt(mineig*maxeig)*((1.0_psb_dpk_/cond + u)/(1.0_psb_dpk_/cond -u))
      c(j) = (zfun(xi(j))/xi(j))*(cn*dn)/((1.0_psb_dpk_/cond - u)**2)
    end do
    eta = -(4.0_psb_dpk_*K*sqrt(mineig*maxeig))/(cond*DPI*N)
    sign = -1.0_psb_dpk_

  end subroutine hhtmethod1

end submodule
