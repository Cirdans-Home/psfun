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
program arnolditest
  ! Test for the parallel computation of matrix function by means of the
  ! psfun_d_arnoldi function. It applies the classical Arnoldi orthogonalization
  ! algorithm on a distributed matrix.
  use psb_base_mod
  use psb_util_mod
  use psfun_d_serial_mod
  use psfun_d_krylov_mod
  implicit none

  ! File input
  character(len=20)           :: mname,fname,filefmt,rhs_file,part,afmt,variant
  real(psb_dpk_)              :: scaling
  ! communications data structure
  type(psb_desc_type):: desc_a
  ! sparse matrices
  type(psb_dspmat_type)  :: a
  type(psb_ldspmat_type) :: aux_a
  ! dense matrix
  real(psb_dpk_), allocatable, target ::  aux_x(:,:)
  real(psb_dpk_), pointer  :: x_col_glob(:)
  type(psb_d_vect_type)    :: x_col, y_col
  ! Matrix function method
  type(psfun_d_serial)        :: fun
  type(psfun_d_krylov)        :: kmethd
  real(psb_dpk_)              :: eps, err
  integer(psb_ipk_)           :: itmax,itrace,istop,iter
  ! blacs parameters
  type(psb_ctxt_type)          :: ctxt
  integer(psb_ipk_)            :: iam, np
  integer(psb_lpk_)            :: lnp
  ! working variables
  integer(psb_ipk_)            :: m_problem, i
  integer(psb_ipk_), allocatable :: ivg(:), perm(:)
  integer(psb_ipk_), allocatable :: ipv(:)
  ! auxiliary parameters
  integer(psb_ipk_)            :: info, ircode
  integer(psb_ipk_), parameter :: iunit=12
  character(len=20)  :: name
  real(psb_dpk_)     :: t1,t2
  integer(psb_epk_)  :: amatsize, descsize, system_size
  real(psb_dpk_), allocatable :: res(:)
  ! real(psb_dpk_), allocatable  :: vy(:)

  info=psb_success_
  name='arnolditest'
  call psb_init(ctxt)
  call psb_info(ctxt,iam,np)
  if (iam < 0) then
    call psb_exit(ctxt) ! This should not happen, but just in case
    stop
  endif
  if(psb_get_errstatus() /= 0) goto 9999
  !
  ! Hello world
  !
  if( iam == psb_root_ ) then
    write(psb_out_unit,'("Welcome to the ",a," program of PSFUN")')name
  end if
  !
  ! Read input information from file
  !
  call get_parms(ctxt,mname,rhs_file,filefmt,part,afmt,fname,&
    & variant,scaling,eps,itmax,itrace,istop)
  if(iam == psb_root_) then
    write(psb_out_unit,*)''
    write(psb_out_unit,'("Solving matrix        : ",a)')mname
    write(psb_out_unit,'("RHS vector            : ",a)')rhs_file
    write(psb_out_unit,'("File format           : ",a)')filefmt
    write(psb_out_unit,'("Partitioning strategy : ",a)')part
    write(psb_out_unit,'("Storage strategy      : ",a)')afmt
    write(psb_out_unit,'("Function to compute   : ",a)')fname
    write(psb_out_unit,'("Algorithmic variant   : ",a)')variant
    write(psb_out_unit,'("Scaling               : ",f8.2)')scaling
  end if

  call psb_barrier(ctxt)

  if (iam == psb_root_) then
    select case(psb_toupper(filefmt))
    case('MM')
      ! For Matrix Market we have an input file for the matrix
      ! and an (optional) second file for the RHS.
      call mm_mat_read(aux_a,info,iunit=iunit,filename=mname)
      if (info == psb_success_) then
        if (rhs_file /= 'NONE') then
          call mm_array_read(aux_x,info,iunit=iunit,filename=rhs_file)
        end if
      end if

    case ('HB')
      ! For Harwell-Boeing we have a single file which may or may not
      ! contain an RHS.
      call hb_read(aux_a,info,iunit=iunit,b=aux_x,filename=mname)

    case default
      info = -1
      write(psb_err_unit,*) 'Wrong choice for fileformat ', filefmt
    end select
    if (info /= psb_success_) then
      write(psb_err_unit,*) 'Error while reading input matrix '
      call psb_abort(ctxt)
    end if

    m_problem = aux_a%get_nrows()
    call psb_bcast(ctxt,m_problem)

    ! At this point aux_x may still be unallocated
    if (size(aux_x,dim=1) == m_problem) then
      ! if any rhs were present, broadcast the first one
      write(psb_err_unit,'("Ok, got an rhs ")')
      x_col_glob =>aux_x(:,1)
      call psb_gelp('N',perm(1:m_problem),&
           & x_col_glob(1:m_problem),info)
    else
      write(psb_out_unit,'("Generating an rhs...")')
      write(psb_out_unit,'(" ")')
      call psb_realloc(m_problem,1,aux_x,ircode)
      if (ircode /= 0) then
         call psb_errpush(psb_err_alloc_dealloc_,name)
         goto 9999
      endif

      x_col_glob => aux_x(:,1)
      do i=1, m_problem
         x_col_glob(i) = done
      enddo
    endif

  else
    call psb_bcast(ctxt,m_problem)

  end if

  ! switch over different partition types
  select case(psb_toupper(part))
  case('BLOCK')
    if (iam == psb_root_) write(psb_out_unit,'("Partition type: block")')
    call psb_matdist(aux_a, a,  ctxt,desc_a,info,fmt=afmt,parts=part_block)

  case('GRAPH')
    if (iam == psb_root_) then
      write(psb_out_unit,'("Partition type: graph vector")')
      write(psb_out_unit,'(" ")')
      !      write(psb_err_unit,'("Build type: graph")')
      call aux_a%cscnv(info,type='csr')
      lnp = np
      call build_mtpart(aux_a,lnp)

    endif
    call psb_barrier(ctxt)
    call distr_mtpart(psb_root_,ctxt)
    call getv_mtpart(ivg)
    call psb_matdist(aux_a, a, ctxt,desc_a,info,fmt=afmt,vg=ivg)

  case default
    if (iam == psb_root_) write(psb_out_unit,'("Partition type: block")')
    call psb_matdist(aux_a, a,  ctxt,desc_a,info,fmt=afmt,parts=part_block)
  end select

  call psb_scatter(x_col_glob,x_col,desc_a,info,root=psb_root_)
  call psb_geall(y_col,desc_a,info)
  call y_col%zero()
  call psb_geasb(y_col,desc_a,info)
  t2 = psb_wtime() - t1


  call psb_amx(ctxt, t2)

  if (iam == psb_root_) then
     write(psb_out_unit,'(" ")')
     write(psb_out_unit,'("Time to read and partition matrix : ",es12.5)')t2
     write(psb_out_unit,'(" ")')
  end if

  call psb_barrier(ctxt)

  ! Set the options for the matrix function
  call fun%set("FNAME",fname,info)
  call fun%set("VARIANT",variant,info)
  call fun%set("SCALING",scaling,info)

  ! Set the options for the Krylov method
  call kmethd%set("KNAME","ARNOLDI",info)

  ! Doing the matrix function computation
  eps = 1e-6
  itmax = 100
  itrace = 1
  istop = 1
  t1 = psb_wtime()
  call kmethd%apply(fun,a,desc_a,y_col,x_col,eps,info,itmax,itrace,istop,iter,err,res)
  t2 = psb_wtime() - t1

  amatsize = a%sizeof()
  descsize = desc_a%sizeof()
  system_size = desc_a%get_global_rows()
  call psb_sum(ctxt,amatsize)
  call psb_sum(ctxt,descsize)

  call psb_barrier(ctxt)

  if (iam == psb_root_) then
    write(psb_out_unit,'(" ")')
    write(psb_out_unit,'("Number of processes           : ",i12)')np
    write(psb_out_unit,'("Matrix size                   : ",i12)')system_size
    write(psb_out_unit,'("Time to solve system          : ",es12.5)')t2
    write(psb_out_unit,'("Time per iteration            : ",es12.5)')t2/iter
    write(psb_out_unit,'("Number of iterations          : ",i12)')iter
    write(psb_out_unit,'("Convergence indicator on exit : ",es12.5)')err
    write(psb_out_unit,'("Info  on exit                 : ",i12)')info
    write(psb_out_unit,'("Total memory occupation for      A: ",i12)')amatsize
    write(psb_out_unit,'("Total memory occupation for DESC_A: ",i12)')descsize
    write(psb_out_unit,'("Storage format for               A: ",a)') a%get_fmt()
    write(psb_out_unit,'("Storage format for          DESC_A: ",a)') desc_a%get_fmt()
  end if

#if defined (WITHGNUPLOTFORTRAN)
  if( iam == psb_root_) &
    & call kmethd%plot(fun,iter,res,info)
#endif

  ! Free the memory
  call psb_gefree(x_col, desc_a, info)
  call psb_gefree(y_col, desc_a, info)
  call psb_spfree(a, desc_a,info)
  call psb_cdfree(desc_a,info)

  call psb_exit(ctxt)
  stop

9999 call psb_error(ctxt)

  stop

  contains

    subroutine  get_parms(ctxt,mname,rhs_file,filefmt,part,afmt,fname,variant,scaling,eps,itmax,itrace,istop)
      ! This subroutine reads the parameters needed to run the serialtest
      ! program from standard input
      type(psb_ctxt_type), intent(in)  :: ctxt
      character(len=*),  intent(out) :: mname,rhs_file,fname,variant,filefmt,part,afmt
      real(psb_dpk_), intent(out)    :: scaling,eps
      integer(psb_ipk_), intent(out) :: itmax, itrace, istop

      integer(psb_ipk_)   :: ip, inp_unit
      integer(psb_ipk_) :: np, iam
      character(len=1024) :: filename

      call psb_info(ctxt, iam, np)

      if (iam == psb_root_) then
        if (command_argument_count()>0) then
          call get_command_argument(1,filename)
          inp_unit = 30
          open(inp_unit,file=filename,action='read',iostat=info)
          if (info /= 0) then
            write(psb_err_unit,*) 'Could not open file ',filename,' for input'
            call psb_abort(ctxt)
            stop
          else
            write(psb_err_unit,*) 'Opened file ',trim(filename),' for input'
          end if
        else
          inp_unit=psb_inp_unit
        end if
        read(inp_unit,*) ip
        if (ip == 8) then
          read(inp_unit,*) mname
          read(inp_unit,*) rhs_file
          read(inp_unit,*) filefmt
          read(inp_unit,*) part
          read(inp_unit,*) afmt
          read(inp_unit,*) fname
          read(inp_unit,*) variant
          read(inp_unit,*) scaling
          eps = 1e-6
          itmax = 100
          itrace = 1
          istop = 1
        else if (ip == 12) then
          read(inp_unit,*) mname
          read(inp_unit,*) rhs_file
          read(inp_unit,*) filefmt
          read(inp_unit,*) part
          read(inp_unit,*) afmt
          read(inp_unit,*) fname
          read(inp_unit,*) variant
          read(inp_unit,*) scaling
          read(inp_unit,*) eps
          read(inp_unit,*) itmax
          read(inp_unit,*) itrace
          read(inp_unit,*) istop
        else
          ! wrong number of parameters on input: all hopes are lost
          call pr_usage(izero)
          call psb_abort(ctxt)
          stop 1
        end if

        if (inp_unit /= psb_inp_unit) then
          close(inp_unit)
        end if
      end if
      ! Broadcast values to all processes
      call psb_bcast(ctxt,mname)
      call psb_bcast(ctxt,rhs_file)
      call psb_bcast(ctxt,filefmt)
      call psb_bcast(ctxt,part)
      call psb_bcast(ctxt,afmt)
      call psb_bcast(ctxt,fname)
      call psb_bcast(ctxt,variant)
      call psb_bcast(ctxt,scaling)
      call psb_bcast(ctxt,eps)
      call psb_bcast(ctxt,itmax)
      call psb_bcast(ctxt,itrace)
      call psb_bcast(ctxt,istop)
      return

    end subroutine get_parms

    subroutine pr_usage(iout)
      ! Prints out information on incorrected program usage
      integer(psb_ipk_) :: iout
      write(iout,*)'incorrect parameter(s) found'
      write(iout,*)' usage:  serialtest mname filefmt fname variant &
           & scaling dump'
      write(iout,*)' where:'
      write(iout,*)'     mname:    Filename of the matrix'
      write(iout,*)'     filefmt:  Format of the matrix (MM,HB)'
      write(iout,*)'     part:     Partitioning strategy (GRAPH,BLOCK)'
      write(iout,*)'     afmt:     Storage format of the matrix (CSR,COO)'
      write(iout,*)'     fname :   Matrix function name'
      write(iout,*)'     variant : Variant of the algorithm for the fname function'
      write(iout,*)'     scaling : Scalar scaling parameter '
    end subroutine pr_usage

end program arnolditest
