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

submodule (psfun_d_krylov_mod) psfun_d_lanczos_mod

contains

  module subroutine psfun_d_lanczos(fun,a,desc_a,y,x,eps,info,itmax,itrace,istop,iter,err)
    use psb_base_mod
    use psfun_d_serial_mod
    implicit none

    type(psfun_d_serial), intent(in)     :: fun  ! Function object
    type(psb_dspmat_type), intent(in)    :: a    ! Distribute sparse matrix
    type(psb_desc_type), intent(in)      :: desc_a ! Descriptor for the sparse matrix
    type(psb_d_vect_type), intent(inout) :: y ! Output vector
    type(psb_d_vect_type), intent(inout) :: x ! Input vector
    real(psb_dpk_), intent(in)           :: eps ! Requested tolerance
    integer(psb_ipk_), intent(out)       :: info  ! Output flag
    integer(psb_ipk_), optional, intent(in)  :: itmax ! Maximum number of iteration
    integer(psb_ipk_), optional, intent(in)  :: itrace ! Trace for logoutput
    integer(psb_ipk_), optional, intent(in)  :: istop ! Stop criterion
    integer(psb_ipk_), optional, intent(out) :: iter ! Number of iteration
    real(psb_dpk_), optional, intent(out) :: err ! Last estimate error

  end subroutine psfun_d_lanczos

end submodule
