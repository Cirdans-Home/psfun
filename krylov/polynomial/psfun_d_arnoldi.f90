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

submodule (psfun_d_krylov_mod) psfun_d_arnoldi_mod

contains

module subroutine psfun_d_arnoldi(fun,a,desc_a,y,x,eps,info,itmax,itrace,istop,iter,err)
  ! Simple polynomial method based on the Arnoldi orthogonalization procedure,
  ! the method builds a basis :math:`V_k` for the Krylov subspace
  !
  ! :math:`\mathcal{K}_k(A,x) = \{x,Ax,\ldots,A^{k-1}x\},`
  !
  ! and approximates :math:`y = f(\alpha A)x \approx \beta_1 V_k f(H_k) e_1`,for
  ! :math:`\beta_1 = \|x\|_2`, :math:`e_1` the first vector of the canonical
  ! base of :math:`\mathbb{R}^k`, and :math:`H_k` the Hessemberg matrix given
  ! by :math:`H_k = V_k^T A V_k`.
  use psb_base_mod
  use psfun_d_serial_mod
  use psb_d_krylov_conv_mod, only: log_header, log_conv, log_end
  implicit none

  type(psfun_d_serial), intent(inout)  :: fun  ! Function object
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


  !Local Variables
  integer(psb_ipk_) :: litmax, itrace_, istop_
  integer(psb_lpk_) :: mglob
  integer(psb_ipk_) :: n_row, n_col, nl, naux, itx, i, k
  real(psb_dpk_)    :: scaling
  !Variables of possible large size:
  real(psb_dpk_), allocatable   :: aux(:)
  real(psb_dpk_), allocatable   :: h(:,:), rs(:), e(:), yk(:)
  type(psb_d_vect_type), allocatable :: v(:)
  type(psb_d_vect_type)              :: w, xt
  ! Error log
  real(psb_dpk_) :: deps, errden, errnum, derr
  integer(psb_ipk_) :: err_act
  integer(psb_ipk_) :: ictxt, iam, np

  character(len=20)           :: name
  character(len=60)           :: methdname

  ! debug
  ! real(psb_dpk_), allocatable    :: va(:)

  ! Information on the distributed environment
  info = psb_success_
  name = 'ARNOLDI'
  methdname = trim(name) // trim("+") // trim(fun%fname) // trim("+") // trim(fun%variant)
  call psb_erractionsave(err_act)
  ictxt = desc_a%get_context()
  call psb_info(ictxt, iam, np)

  mglob = desc_a%get_global_rows()
  n_row = desc_a%get_local_rows()
  n_col = desc_a%get_local_cols()

  ! Select the stop-criterion
  if (present(istop)) then
    istop_ = istop
  else
    istop_ = 1
  endif

  if ((istop_ < 1 ).or.(istop_ > 2 ) ) then
    info=psb_err_invalid_istop_
    err=info
    call psb_errpush(info,name,i_err=(/istop_/))
    goto 9999
  endif

  ! Maximum number of iterations
  if (present(itmax)) then
    litmax = itmax
  else
    litmax = 200
  endif
  nl = litmax

  naux = 4*n_col

  ! What do i print?
  if (present(itrace)) then
    itrace_ = itrace
  else
    itrace_ = 0
  end if

  ! Chek sanity of the inputs both on the local and the global
  if (.not.allocated(x%v)) then
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  if (.not.allocated(y%v)) then
    info = psb_err_invalid_vect_state_
    call psb_errpush(info,name)
    goto 9999
  endif
  call psb_chkvect(mglob,lone,x%get_nrows(),lone,lone,desc_a,info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='psb_chkvect on X')
    goto 9999
  end if
  call psb_chkvect(mglob,lone,y%get_nrows(),lone,lone,desc_a,info)
  if(info /= psb_success_) then
    info=psb_err_from_subroutine_
    call psb_errpush(info,name,a_err='psb_chkvect on Y')
    goto 9999
  end if

  ! Allocate memory for the Krylov basis and the auxiliary quantities
  allocate(aux(naux),h(nl+1,nl+1),rs(nl+1),stat=info)

  if (info == psb_success_) call psb_geall(v,desc_a,info,n=nl+1)
  if (info == psb_success_) call psb_geall(w,desc_a,info)
  if (info == psb_success_) call psb_geall(xt,desc_a,info)
  if (info == psb_success_) call psb_geasb(v,desc_a,info,mold=x%v)
  if (info == psb_success_) call psb_geasb(w,desc_a,info,mold=x%v)
  if (info == psb_success_) call psb_geasb(xt,desc_a,info,mold=x%v)
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_non_
    call psb_errpush(info,name)
    goto 9999
  end if

  ! We use the scaling in the computation of the Krylov subspace
  scaling = fun%scaling
  call fun%set("SCALING",done,info)

  ! LogVariables
  errnum = dzero
  errden = done
  deps   = eps
  ! Start the iteration
  if ((itrace_ > 0).and.(iam == 0)) call log_header(methdname)

  call psb_geaxpby(done,x,dzero,v(1),desc_a,info) ! v(1) = x
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_non_
    call psb_errpush(info,name)
    goto 9999
  end if
  rs(1) = psb_genrm2(v(1),desc_a,info) ! rs(1) = ||v(1)||_2
  rs(2:) = dzero
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_non_
    call psb_errpush(info,name)
    goto 9999
  end if

  call v(1)%scal(done/rs(1)) !v(1) = v(1)/||v(1)||_2
  itx = 0
  arnoldicycle: do i = 2, nl, 1
    itx = itx + 1
    ! write(*,'("Iteration ",i2)')i
    call psb_spmm(scaling,a,v(i-1),dzero,w,desc_a,info,work=aux) ! w = A v(i)
    ! va = w%get_vect()
    ! write(*,'("w(1) =",f8.2)')va(1)
    if (info /= psb_success_) then
      info=psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
    end if
    do k = 1, i-1, 1
      h(k,i-1) = psb_gedot(v(k),w,desc_a,info)            !h_{k,i-1}=w^T*v(k)
      ! write(*,'("h(",i2,","i2") =",f8.2)')k,i-1,h(k,i-1)
      call psb_geaxpby(-h(k,i-1),v(k),done,w,desc_a,info) !w = w - h_{k,i-1}v(k)
    end do
    h(i,i-1) = psb_genrm2(w,desc_a,info)                  !h_{i,i-1} = ||w||_2
    ! write(*,'("h(",i2,","i2") =",f8.2)')i,i-1,h(i,i-1)
    call psb_geaxpby(done/h(i,i-1),w,dzero,v(i),desc_a,info) !v(i) = w/h_{i,i-1}

    ! A posteriori error estimate:
    ! This is the default a posteriori-error estimate, it costs the computation
    ! of the matrix function at each iteration to use the computed values
    if( istop_ == 1) then
      !
      ! build y and then compute the residual as
      ! rs(j) = \| x \|_2 h_{i,i-1} | e_{i-1}^T f(H_{i-1}) e_1 |
      !
      allocate(e(i-1), stat=info)
      if (info /= 0) then
        info=psb_err_from_subroutine_non_
        call psb_errpush(info,name)
        goto 9999
      end if
      e(:) = dzero

      allocate(yk(i-1), stat=info)
      if (info /= 0) then
        info=psb_err_from_subroutine_non_
        call psb_errpush(info,name)
        goto 9999
      end if

      e(1) = done
      call fun%apply(h(1:i-1,1:i-1),yk,e,info)
      if (info /= psb_success_) then
        info=psb_err_from_subroutine_non_
        call psb_errpush(info,name)
        goto 9999
      end if

      errnum = rs(1)*h(i,i-1)*abs( yk(1) );
      errden = norm2(yk)
      rs(i) = errnum/errden

      ! log of the error estimate
      if (itrace_ > 0) &
           & call log_conv(methdname,iam,itx,itrace_,errnum,errden,deps)

      if (allocated(e)) deallocate(e, stat=info)
      if (info /= 0) then
        info=psb_err_from_subroutine_non_
        call psb_errpush(info,name)
        goto 9999
      end if
      ! Check convergence
      if(rs(i) < deps) then
        exit arnoldicycle
      else
        if (allocated(yk)) deallocate(yk, stat=info)
        if (info /= 0) then
          info=psb_err_from_subroutine_non_
          call psb_errpush(info,name)
          goto 9999
        end if
      end if

    end if

  end do arnoldicycle
  !
  ! Assemble the final solution
  !
  if(.not.allocated(yk)) then
    ! First we look if we have arrived here having not used an a posteriori
    ! error estimate, and thus not having computed the solution in the Krylov
    ! subspace.
    allocate(e(itx-1), stat=info)
    if (info /= 0) then
      info=psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
    end if
    e(:) = dzero

    allocate(yk(itx-1), stat=info)
    if (info /= 0) then
      info=psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
    end if
    e(1) = done
    call fun%apply(h(1:itx-1,1:itx-1),yk,e,info)
    if (info /= psb_success_) then
      info=psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
    end if

    if (allocated(e)) deallocate(e, stat=info)
    if (info /= 0) then
      info=psb_err_from_subroutine_non_
      call psb_errpush(info,name)
      goto 9999
    end if
  end if

  ! Final Log
  call log_end(methdname,iam,itx,itrace_,errnum,errden,deps,err=derr,iter=iter)
  if (present(err)) err = derr

  call  psb_geaxpby(yk(1),v(1),dzero,y,desc_a,info)
  do i=2,itx-1,1
    call  psb_geaxpby(yk(i),v(i),done,y,desc_a,info)
  end do
  call y%scal(rs(1))

  ! Set back the scaling
  call fun%set("SCALING",scaling,info)

  ! Free the memory
  if (info == psb_success_) call psb_gefree(v,desc_a,info)
  if (info == psb_success_) call psb_gefree(w,desc_a,info)
  if (info == psb_success_) call psb_gefree(xt,desc_a,info)
  if (info == psb_success_) deallocate(aux,h,rs,stat=info)
  if (info /= psb_success_) then
    info=psb_err_from_subroutine_non_
    call psb_errpush(info,name)
    goto 9999
  end if

  ! Close and return
  call psb_erractionrestore(err_act)
return

9999 call psb_error_handler(err_act)
return

end subroutine psfun_d_arnoldi

end submodule
