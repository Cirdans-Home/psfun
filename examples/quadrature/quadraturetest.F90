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
program quadraturetest
  use psb_base_mod
  use psfun_utils_mod
  use psfun_quadrature_mod

  type(psb_ctxt_type)      :: ctxt
  type(psfun_z_quadrature) :: quad
  procedure (zquadfun), pointer  :: zfun
  procedure (zquadrule), pointer :: quadformula
  integer(psb_ipk_)        :: N, info
  real(psb_dpk_)           :: rparams(2)
  character(len=20)        :: name
  ! Variable for debug
  complex(psb_dpk_), allocatable, dimension(:) :: xi     ! Poles of the formula
  complex(psb_dpk_), allocatable, dimension(:) :: c      ! Scaling of the formula
  real(psb_dpk_)    :: eta! Global Scaling
  real(psb_dpk_)    :: sign   ! Sign for A
  integer(psb_ipk_) :: i

  info=psb_success_
  name='quadraturetest'
  call psb_init(ctxt)
  call psb_info(ctxt,iam,np)
  if (iam < 0) then
    call psb_exit(ctxt) ! This should not happen, but just in case
    stop
  endif
  if(psb_get_errstatus() /= 0) goto 9999

  N = 10
  rparams(1) = 0.01_psb_dpk_
  rparams(2) = 4.0_psb_dpk_

  zfun => fun
  quadformula => hhtmethod1
  call quad%computepoles(quadformula=quadformula,&
    & zfun=zfun,N=N,info=info,rparams=rparams)
  call quad%plot(zfun,info)

  call psb_exit(ctxt)
  stop

9999 call psb_error(ctxt)

  stop

contains

  function fun(z) result(res)
    use psb_base_mod
    implicit none
    complex(psb_dpk_), intent(in) :: z
    complex(psb_dpk_) :: res

    res = sqrt(1.0 + z)

  end function

end program
