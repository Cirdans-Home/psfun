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
    !! The [[psfun_d_krylov_mod]] contains the generic call to a Krylov subspace
    !! method for the computation of \(y = f(A) x\), for \(A\) large and
    !! sparse.

  use psb_base_mod
  use psb_krylov_mod
  use amg_prec_mod
  use psfun_d_serial_mod
  use ogpf
  implicit none

  type, public :: psfun_d_krylov
    character(len=20)     :: kname   = 'ARNOLDI' !! Name of the Krylov method
    character(len=20)     :: kmethd !! Method for the solution of the linear system
    type(amg_dprec_type)  :: prec   !! Preconditioner for the inner solution method
  contains
    ! Set the options
    procedure, pass(meth) :: setstring  => psfun_d_setstring
    generic, public :: set => setstring
    procedure, pass(meth) :: apply      => psfun_d_parallel_apply
    procedure, pass(meth) :: plot       => psfun_d_plot_info
    procedure, pass(meth) :: precinit   => psfun_d_prec_init
    procedure, pass(meth) :: precbuild  => psfun_d_prec_build
  end type

  private :: psfun_d_setstring

    !! The various Krylov methods are contained in associated submodules, the
    !! idea is to have different methods, associated to the same
    !! orthogonalization method in the same submodule
  interface
      !! Simple polynomial method based on the Arnoldi orthogonalization procedure,
      !! the method builds a basis \(V_k\) for the Krylov subspace
      module subroutine psfun_d_arnoldi(fun,a,desc_a,y,x,eps,info,itmax,itrace,istop,iter,err,res)
          !! Simple polynomial method based on the Arnoldi orthogonalization procedure,
          !! the method builds a basis \(V_k\) for the Krylov subspace
          !
          !! \begin{equation*}
          !! \mathcal{K}_k(A,x) = \{x,Ax,\ldots,A^{k-1}x\},
          !! \end{equation*}
          !
          !! and approximates \(y = f(\alpha A)x \approx \beta_1 V_k f(H_k) e_1\),for
          !! \(\beta_1 = \|x\|_2\), \(e_1\) the first vector of the canonical
          !! base of \(\mathbb{R}^k\), and \(H_k\) the Hessemberg matrix given
          !! by \(H_k = V_k^T A V_k\).
          type(psfun_d_serial), intent(inout)  :: fun  !! Function object
          type(psb_dspmat_type), intent(in)    :: a    !! Distribute sparse matrix
          type(psb_desc_type), intent(in)      :: desc_a !! Descriptor for the sparse matrix
          type(psb_d_vect_type), intent(inout) :: y !! Output vector
          type(psb_d_vect_type), intent(inout) :: x !! Input vector
          real(psb_dpk_), intent(in)           :: eps !! Requested tolerance
          integer(psb_ipk_), intent(out)       :: info  !! Output flag
          integer(psb_ipk_), optional, intent(in)  :: itmax !! Maximum number of iteration
          integer(psb_ipk_), optional, intent(in)  :: itrace !! Trace for logoutput
          integer(psb_ipk_), optional, intent(in)  :: istop !! Stop criterion
          integer(psb_ipk_), optional, intent(out) :: iter !! Number of iteration
          real(psb_dpk_), optional, intent(out) :: err !! Last estimate error
          real(psb_dpk_), optional, allocatable, intent(out) :: res(:) !! Vector of the residuals
      end subroutine
  end interface

  interface
      !! Simple polynomial method based on the Lanczos orthogonalization procedure,
      !! the method builds a basis \(V_k\) for the Krylov subspace
      module subroutine psfun_d_lanczos(fun,a,desc_a,y,x,eps,info,itmax,itrace,istop,iter,err,res)
        !! Simple polynomial method based on the Lanczos orthogonalization procedure,
        !! the method builds a basis \(V_k\) for the Krylov subspace
        !!
        !! \begin{equation*}
        !! \mathcal{K}_k(A,x) = \{x,Ax,\ldots,A^{k-1}x\},
        !! \end{equation*}
        !!
        !! and approximates \(y = f(\alpha A)x \approx \beta_1 V_k f(T_k) e_1\),for
        !! \(\beta_1 = \|x\|_2\), \(e_1\) the first vector of the canonical
        !! base of \(\mathbb{R}^k\), and \(T_k\) the Symmetric tridiagonal
        !! matrix given by \(T_k = V_k^T A V_k\).
        type(psfun_d_serial), intent(inout)  :: fun  !! Function object
        type(psb_dspmat_type), intent(in)    :: a    !! Distribute sparse matrix
        type(psb_desc_type), intent(in)      :: desc_a !! Descriptor for the sparse matrix
        type(psb_d_vect_type), intent(inout) :: y !! Output vector
        type(psb_d_vect_type), intent(inout) :: x !! Input vector
        real(psb_dpk_), intent(in)           :: eps !! Requested tolerance
        integer(psb_ipk_), intent(out)       :: info  !! Output flag
        integer(psb_ipk_), optional, intent(in)  :: itmax !! Maximum number of iteration
        integer(psb_ipk_), optional, intent(in)  :: itrace !! Trace for logoutput
        integer(psb_ipk_), optional, intent(in)  :: istop !! Stop criterion
        integer(psb_ipk_), optional, intent(out) :: iter !! Number of iteration
        real(psb_dpk_), optional, intent(out) :: err !! Last estimate error
        real(psb_dpk_), optional, allocatable, intent(out) :: res(:) !! Vector of the residuals
      end subroutine psfun_d_lanczos
  end interface

contains

  subroutine psfun_d_setstring(meth,what,val,info)
      !! Set function for setting options defined by a string
      use psb_base_mod
      implicit none

      class(psfun_d_krylov), intent(inout) :: meth  ! Krylov method object
      character(len=*), intent(in)         :: what  ! String of option to set
      character(len=*), intent(in)         :: val   ! Value of the string
      integer(psb_ipk_), intent(out)       :: info  ! Output flag

      info = psb_success_
      select case (psb_toupper(what))
      case ("KNAME")
          meth%kname = val
      case default
          info = psb_err_invalid_args_combination_
      end select

  end subroutine psfun_d_setstring

  subroutine psfun_d_parallel_apply(meth,fun,a,desc_a,y,x,eps,info,itmax,itrace,istop,iter,err,res)
      !! This is the generic function for applying every implemented Krylov
      !! method. The general iteration parameters (like the number of iteration,
      !! the stop criterion to be used, and the verbosity of the trace) can be
      !! passed directly to this routine. All the constitutive parameters of
      !! the actual method, and the information relative to the function are
      !! instead contained in the meth and fun objects. The Descriptor object
      !! `desc_a' contains the properties of the parallel environment.
      use psb_base_mod
      use psfun_d_serial_mod
      implicit none

      class(psfun_d_krylov), intent(in)    :: meth !! Krylov method object
      type(psfun_d_serial), intent(inout)  :: fun  !! Function object
      type(psb_dspmat_type), intent(in)    :: a    !! Distribute sparse matrix
      type(psb_desc_type), intent(in)      :: desc_a !! Descriptor for the sparse matrix
      type(psb_d_vect_type), intent(inout) :: y !! Output vector
      type(psb_d_vect_type), intent(inout) :: x !! Input vector
      real(psb_dpk_), intent(in)           :: eps !! Requested tolerance
      integer(psb_ipk_), intent(out)       :: info  !! Output flag
      integer(psb_ipk_), optional, intent(in)  :: itmax !! Maximum number of iteration
      integer(psb_ipk_), optional, intent(in)  :: itrace !! Trace for logoutput
      integer(psb_ipk_), optional, intent(in)  :: istop !! Stop criterion
      integer(psb_ipk_), optional, intent(out) :: iter !! Number of iteration
      real(psb_dpk_), optional, intent(out) :: err !! Last estimate error
      real(psb_dpk_), optional, allocatable, intent(out) :: res(:) !! Vector of the residuals

      select case (psb_toupper(meth%kname))
      case ("ARNOLDI")
          call psfun_d_arnoldi(fun,a,desc_a,y,x,eps,info,itmax,itrace,istop,iter,err,res)
      case ("LANCZOS")
          call psfun_d_lanczos(fun,a,desc_a,y,x,eps,info,itmax,itrace,istop,iter,err,res)
      case default
          info = psb_err_invalid_args_combination_
      end select

  end subroutine psfun_d_parallel_apply

  subroutine psfun_d_plot_info(meth,fun,iter,res,info)
      !! This function plots the convergence history of the Krylov method
      use psb_base_mod
      use ogpf

      implicit none

      class(psfun_d_krylov), intent(inout)      :: meth !! Krylov method
      type(psfun_d_serial), intent(inout)       :: fun  !! Function object
      integer(psb_ipk_), intent(in)             :: iter !! Number of iteration
      real(psb_dpk_), intent(in), dimension(:)  :: res  !! Residual vector
      integer(psb_ipk_), intent(out)            :: info !! Result of the Gnuplot call

      type(gpf) :: gp
      character(len=100) :: pf
      real(psb_dpk_), allocatable :: xarray(:)
      integer(psb_ipk_) :: n, indplot
      real(psb_dpk_) :: plot

      n = size(res)
      if (n <= iter+1) then
          indplot = n
          plot = n
      else
          indplot = iter+1
          plot = iter+1
      end if

#if defined(WITHGNUPLOTFORTRAN)
      allocate(xarray(indplot), stat=info)
      xarray = linspace(1.0_psb_dpk_,plot,indplot)

      pf = trim(meth%kname) // trim("-") // trim(fun%fname) &
            & // trim("-") // trim(fun%variant) // trim("-Residual")
      call gp%options("set terminal pdf;&
                    &set output '"//trim(pf)//".pdf'")
      call gp%title("Residual history "// trim(meth%kname) // trim("-") &
        & // trim(fun%fname) // trim("-") // trim(fun%variant))
      call gp%xlabel('Iteration')
      call gp%ylabel('Relative Residual')
      call gp%options('set style data linespoints;&
                    &set logscale y;&
                    &set format y "10^{%L}";')
      !Call Plot to draw a vector against a vector of data
      call gp%plot(xarray(1:indplot), res(1:indplot),'with lp lt 7')
      info = psb_success_

      if (allocated(xarray)) deallocate(xarray, stat=info)
#else
      info = psb_err_from_subroutine_
#endif


  end subroutine psfun_d_plot_info

  subroutine psfun_d_prec_init(meth,ctxt,ptype,info)
      !! This function performs the init of the preconditioner for the inner
      !! solve in the rational Krylov method
      use psb_base_mod
      use amg_prec_mod

      implicit none

      class(psfun_d_krylov), intent(inout)      :: meth  !! Krylov method
      type(psb_ctxt_type), intent(in)           :: ctxt  !! Parallel context
      character(len=20), intent(in)             :: ptype !! PSBLAS/AMG4PSBLAS preconditioner
      integer(psb_ipk_)                         :: info  !! Result of the init call

      ! Local Variables
      integer(psb_ipk_) :: err_act
      character(len=20) :: name

      name = 'd_prec_init'
      call psb_erractionsave(err_act)

      info = psb_success_

      call meth%prec%init(ctxt,ptype,info)

      if(info /= psb_success_) then
          call psb_errpush(info,name)
          goto 9999
      end if

      call psb_erractionrestore(err_act)
      return

9999 call psb_error_handler(err_act)
    return

  end subroutine psfun_d_prec_init

  subroutine psfun_d_prec_build(meth,a,desc_a,info)
      !! This function builds the AMG4PSBLAS preconditioner for the inner solve
      !! in a Rational Krylov method
      use psb_base_mod
      use amg_prec_mod

      implicit none

      class(psfun_d_krylov), intent(inout)      :: meth   !! Krylov method
      type(psb_dspmat_type), intent(inout)      :: a      !! Sparse matrix
      type(psb_desc_type), intent(inout)        :: desc_a !! Descriptor for the sparse matrix
      integer(psb_ipk_)                         :: info   !! Result of the init call

      ! Local Variables
      integer(psb_ipk_) :: err_act
      character(len=20) :: name

      name = 'd_prec_build'
      call psb_erractionsave(err_act)

      ! Build The Hierarchy
      call meth%prec%hierarchy_build(a,desc_a,info)
      if(info /= psb_success_) then
          call psb_errpush(info,name)
          goto 9999
      end if

      ! Build The Smoothers
      call meth%prec%smoothers_build(a,desc_a,info)
      if(info /= psb_success_) then
          call psb_errpush(info,name)
          goto 9999
      end if

      call psb_erractionrestore(err_act)
      return

9999 call psb_error_handler(err_act)
    return

  end subroutine psfun_d_prec_build

end module psfun_d_krylov_mod
