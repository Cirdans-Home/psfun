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

module psfun_krylov_mod
  !! This module implements the Krylov methods for the solution of shifted linear
  !! systems: \((η A + ζ I)x = b\). It is built as an extension of the Krylov module
  !! for PSBLAS (`psb_krylov_mod`)
  use psb_base_mod
  use psb_prec_mod

interface psb_krylov
  !! We add the methods for the shifted system to the same interfaces in PSBLAS
  !! it is a lazy way to avoid modifying all the Krylov methods in PSBLAS to
  !! allow for the solution of shifted linear systems
  module procedure psfun_dkrylov_vect
end interface psb_krylov

contains

  subroutine psfun_dkrylov_vect(method,a,prec,b,eta,zeta,x,eps,desc_a,info,&
       & itmax,iter,err,itrace,irst,istop,cond)
    !! Apply Krylov method to \( (\eta A + \zeta I)x = b\) on distributed vectors
    use psb_base_mod, only  : psb_ipk_, psb_desc_type, psb_dspmat_type, &
         & psb_dpk_, psb_d_vect_type
    use psb_prec_mod, only : psb_dprec_type

    character(len=*)                      :: method
    Type(psb_dspmat_type), Intent(in)     :: a
    Type(psb_desc_type), Intent(in)       :: desc_a
    class(psb_dprec_type), intent(inout)  :: prec
    type(psb_d_vect_type), Intent(inout)  :: b
    real(psb_dpk_), intent(in)            :: eta
    real(psb_dpk_), intent(in)            :: zeta
    type(psb_d_vect_type), Intent(inout)  :: x
    Real(psb_dpk_), Intent(in)            :: eps
    integer(psb_ipk_), intent(out)                  :: info
    integer(psb_ipk_), Optional, Intent(in)         :: itmax, itrace, irst,istop
    integer(psb_ipk_), Optional, Intent(out)        :: iter
    Real(psb_dpk_), Optional, Intent(out) :: err,cond

    abstract interface
      subroutine psfun_dkryl_cond_vect(a,prec,b,eta,zeta,x,eps,desc_a,info,&
           &itmax,iter,err, itrace,istop,cond)
        import :: psb_ipk_, psb_dpk_, psb_desc_type, &
             & psb_dspmat_type, psb_dprec_type, psb_d_vect_type
        Type(psb_dspmat_type), Intent(in)    :: a
        Type(psb_desc_type), Intent(in)      :: desc_a
        class(psb_dprec_type), intent(inout) :: prec
        type(psb_d_vect_type), Intent(inout) :: b
        real(psb_dpk_), intent(in)            :: eta
        real(psb_dpk_), intent(in)            :: zeta
        type(psb_d_vect_type), Intent(inout) :: x
        Real(psb_dpk_), Intent(in)           :: eps
        integer(psb_ipk_), intent(out)                 :: info
        integer(psb_ipk_), Optional, Intent(in)        :: itmax, itrace,istop
        integer(psb_ipk_), Optional, Intent(out)       :: iter
        Real(psb_dpk_), Optional, Intent(out) :: err, cond
      end subroutine psfun_dkryl_cond_vect
    end interface

    procedure(psfun_dkryl_cond_vect) :: psfun_dcg_vect

    logical           :: do_alloc_wrk
    type(psb_ctxt_type) :: ctxt
    integer(psb_ipk_) :: me,np,err_act, itrace_
    character(len=20)             :: name

    info = psb_success_
    name = 'psfun_krylov'
    call psb_erractionsave(err_act)

    ctxt=desc_a%get_context()

    call psb_info(ctxt, me, np)

    ! Default return for COND
    if (present(cond)) cond = dzero

    if (present(itrace)) then
      itrace_ = itrace
    else
      itrace_ = -1
    end if

    do_alloc_wrk = .not.prec%is_allocated_wrk()
    if (do_alloc_wrk) call prec%allocate_wrk(info,vmold=x%v,desc=desc_a)

    select case(psb_toupper(method))
    case('CG')
      call  psfun_dcg_vect(a,prec,b,eta,zeta,x,eps,desc_a,info,&
           &itmax,iter,err,itrace=itrace_,istop=istop,cond=cond)
      case default
       if (me == 0) write(psb_err_unit,*) trim(name),&
            & ': Warning: Unknown method  ',method,&
            & ', defaulting to CG'
       call  psfun_dcg_vect(a,prec,b,eta,zeta,x,eps,desc_a,info,&
            &itmax,iter,err,itrace=itrace_,istop=istop)
    end select

    if ((info==psb_success_).and.do_alloc_wrk) call prec%free_wrk(info)

    if(info /= psb_success_) then
     info = psb_err_from_subroutine_
     call psb_errpush(info,name,a_err=trim(method))
     goto 9999
    end if

    call psb_erractionrestore(err_act)
    return

9999 call psb_error_handler(ctxt,err_act)
    return

  end subroutine psfun_dkrylov_vect

end module psfun_krylov_mod
