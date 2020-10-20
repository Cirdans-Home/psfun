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
program utiltest

  use psb_base_mod
  use psfun_utils_mod
  implicit none

  ! blacs parameters
  integer(psb_ipk_)           :: ictxt, iam, np
  ! flags
  integer(psb_ipk_)           :: info
  ! Variables
  real(psb_dpk_)              :: K(2),L,u,sn,cn,dn,result
  real(psb_dpk_), parameter :: coeffs(7) = (/132.0,42.0,14.0,5.0,2.0,1.0,0.0/)
  complex(psb_dpk_)           :: cu,csn,ccn,cdn

  info=psb_success_
  call psb_init(ictxt)
  call psb_info(ictxt,iam,np)

  if (iam < 0) then
    call psb_exit(ictxt) ! This should not happen, but just in case
    stop
  endif
  if(psb_get_errstatus() /= 0) goto 9999

  if( (iam == 0).and.(np > 1) ) then
    write(psb_err_unit, *) "This is a serial example, number of processes is ",np
    goto 9999
  else if( (iam == 0).and.(np == 1)) then
    write(psb_out_unit,*) "Welcome to the utilstest program of PSFUN"
  end if

  ! Test of Horner rule for polynomial evaluation
  write(psb_out_unit,*)
  write(psb_out_unit,'("Polynomial Evaluation")')
  result = horner(coeffs,1.0_psb_dpk_)
  write(psb_out_unit,'("p(1) = ",f17.0," (should be 196)")')result

  ! Test of Elliptic Jacobi Integrals and Functions
  write(psb_out_unit,*)
  write(psb_out_unit,'("Elliptic Integrals and Jacobi Functions")')
  L = 0.5_psb_dpk_
  K = ellipkkp(L)
  write(psb_out_unit,'("L = ",f3.1," K = ",f17.15," Kp = ",f17.15)')L,K(1),K(2)
  write(psb_out_unit,'("L = 0.5 K = 1.588191701877384 Kp = 2.978718299395645 (Control Values)")')
  write(psb_out_unit,'("  Jacobi Elliptic Functions")')
  u = 1.0_psb_dpk_
  L = 0.5_psb_dpk_
  call ellipj(u,L,sn,cn,dn)
  write(psb_out_unit,'("sn = ",f17.15," cn = ",f17.15," dn = ",f17.15)')sn,cn,dn
  write(psb_out_unit,'("sn = 0.838274911024583 cn = 0.545247809300255 dn = 0.984699634947678 (Control Values)")')
  write(psb_out_unit,'("  Jacobi Elliptic Functions (Complex)")')
  cu = cmplx(0.0,0.5_psb_dpk_)
  L = 1.0_psb_dpk_
  call ellipj(cu,L,csn,ccn,cdn)
  write(psb_out_unit,'("sn = (",f3.1,",",f17.15,") cn = (",f17.15,",",f3.1,") dn = (",f17.15,",",f3.1,")")')csn,ccn,cdn
  write(psb_out_unit,'("sn = (0.0,0.521141424070547) cn = (1.127647278133671,0.0) dn = (1.000253555731494,0.0) (Control Values)")')
  call psb_exit(ictxt)
  stop

9999 call psb_error(ictxt)

  stop

end program utiltest
