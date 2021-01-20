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
    ! The psfun_d_krylov_mod contains the generic call to a Krylov subspace
    ! method for the computation of :math:`y = f(A) x`, for :math:`A` large and
    ! sparse.

  use psb_base_mod
  use psfun_d_serial_mod
#if defined(WITHGNUPLOTFORTRAN)
  ! If the library is compiled with the GNUPLOT Fortran library option we
  ! include here the relative modules that can be used to plot the convergence
  ! history of the Krylov method
  use datatypes
  use gnuplot_module_data
  use gnuplot_module
#endif
  implicit none

  type, public :: psfun_d_krylov
    character(len=20)   :: kname   = 'ARNOLDI' ! Name of the Krylov method
  contains
    ! Set the options
    procedure, pass(meth) :: setstring  => psfun_d_setstring
    generic, public :: set => setstring
    procedure, pass(meth) :: apply      => psfun_d_parallel_apply
#if defined(WITHGNUPLOTFORTRAN)
    procedure, pass(meth) :: plot       => psfun_d_plot_info
#endif
  end type

  private :: psfun_d_setstring

    ! The various Krylov methods are contained in associated submodules, the
    ! idea is to have different methods, associated to the same
    ! orthogonalization method in the same submodule
  interface
      module subroutine psfun_d_arnoldi(fun,a,desc_a,y,x,eps,info,itmax,itrace,istop,iter,err,res)
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
          real(psb_dpk_), optional, allocatable, intent(out) :: res(:) ! Vector of the residuals
      end subroutine
  end interface

  interface
      module subroutine psfun_d_lanczos(fun,a,desc_a,y,x,eps,info,itmax,itrace,istop,iter,err,res)
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
        real(psb_dpk_), optional, allocatable, intent(out) :: res(:) ! Vector of the residuals
      end subroutine psfun_d_lanczos
  end interface

contains

  subroutine psfun_d_setstring(meth,what,val,info)
      ! Set function for setting options defined by a string
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
      ! This is the generic function for applying every implemented Krylov
      ! method. The general iteration parameters (like the number of iteration,
      ! the stop criterion to be used, and the verbosity of the trace) can be
      ! passed directly to this routine. All the constitutive parameters of
      ! the actual method, and the information relative to the function are
      ! instead contained in the meth and fun objects. The Descriptor object
      ! desc_a contains the properties of the parallel environment.
      use psb_base_mod
      use psfun_d_serial_mod
      implicit none

      class(psfun_d_krylov), intent(in)    :: meth ! Krylov method object
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
      real(psb_dpk_), optional, allocatable, intent(out) :: res(:) ! Vector of the residuals

      select case (psb_toupper(meth%kname))
      case ("ARNOLDI")
          call psfun_d_arnoldi(fun,a,desc_a,y,x,eps,info,itmax,itrace,istop,iter,err,res)
      case ("LANCZOS")
          call psfun_d_lanczos(fun,a,desc_a,y,x,eps,info,itmax,itrace,istop,iter,err,res)
      case default
          info = psb_err_invalid_args_combination_
      end select

  end subroutine psfun_d_parallel_apply

#if defined(WITHGNUPLOTFORTRAN)
  subroutine psfun_d_plot_info(meth,fun,iter,res,info)
      ! This function plots the convergence history of the Krylov method
      use psb_base_mod
      use psfun_utils_mod, only: linspace
      use datatypes
      use gnuplot_module_data
      use gnuplot_module

      implicit none

      class(psfun_d_krylov), intent(inout)      :: meth ! Krylov method
      type(psfun_d_serial), intent(inout)       :: fun  ! Function object
      integer(psb_ipk_), intent(in)             :: iter ! Number of iteration
      real(psb_dpk_), intent(in), dimension(:)  :: res  ! Residual vector
      integer(psb_ipk_), intent(out)            :: info ! Result of the Gnuplot call

      ! local variable
      type(gnuplot_ctrl), pointer :: ptr_gctrl
      character(len=GP_CMD_SIZE)  :: cmd
      character(len=100) :: pf
      real(psb_dpk_), allocatable :: xarray(:)
      integer(psb_ipk_) :: n, indplot
      real(psb_dpk_) :: plotn

      GNUPLOT_SHOWWARNINGS=.false.
      n = size(res)
      if (n <= iter) then
          plotn = n
          indplot = n
      else
          plotn = iter
          indplot = iter
      end if

      allocate(xarray(indplot), stat=info)
      call linspace(1.0_psb_dpk_,plotn,xarray )

      pf = trim(meth%kname) // trim("-") // trim(fun%fname) &
            & // trim("-") // trim(fun%variant) // trim("-Residual")

      ! Actual plot
      ptr_gctrl=>gnuplot_init('-clear')
      info = gnuplot_hardcopy(ptr_gctrl,'PNG',trim(pf)//trim(".png"))
      info = gnuplot_setrange(ptr_gctrl,'x',1.0_psb_dpk_,plotn)
      info = gnuplot_setrange(ptr_gctrl,'y',minval(res(1:indplot)),maxval(res(1:indplot)))
      info = gnuplot_settitle(ptr_gctrl,pf)
      info = gnuplot_setaxislabel(ptr_gctrl,'x','Iteration')
      info = gnuplot_setaxislabel(ptr_gctrl,'y','Residual')
      info = gnuplot_setscale(ptr_gctrl,'x','NLG')
      info = gnuplot_setscale(ptr_gctrl,'y','LOG')
      info = gnuplot_setstyle(ptr_gctrl,'linespoints')
      info = gnuplot_plot2d(ptr_gctrl,indplot,xarray(1:indplot),res(1:indplot),'residual')
      cmd  = 'set output'
      info = gnuplot_cmd(ptr_gctrl,cmd)
      info = gnuplot_close(ptr_gctrl)

     if (allocated(xarray)) deallocate(xarray, stat=info)

  end subroutine psfun_d_plot_info
#endif

end module psfun_d_krylov_mod
