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

program serialtest
  ! Test program for the serial part of the library. This test program loads a
  ! matrix from file together with some options to test the serial computation
  ! of the matrix functions. Substantially, it test the interfacing with the
  ! library doing the serial part.
  use psb_base_mod
  use psfun_d_serial_mod
  use psb_util_mod, only: mm_mat_read, mm_array_write
  implicit none

  ! File input
  character(len=20)           :: mname,fname,variant
  real(psb_dpk_)              :: scaling
  logical                     :: dump
  ! blacs parameters
  integer(psb_ipk_)           :: ictxt, iam, np
  ! Matrices and vectors
  type(psb_dspmat_type)       :: a ! Needed for reading from sparse matrix
  real(psb_dpk_), allocatable :: x(:),y(:)
  integer(psb_ipk_)           :: n, nnz
  ! Matrix function method
  type(psfun_d_serial)        :: fun
  ! auxiliary parameters
  integer(psb_ipk_)            :: info
  integer(psb_ipk_), parameter :: iunit=12
  real(psb_dpk_)               :: t1,t2

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
    write(psb_out_unit,*) "Welcome to the serialtest program of PSFUN"
  end if

  ! Read input information from file
  call get_parms(ictxt,mname,fname,variant,scaling,dump)
  write(psb_out_unit,*)''
  write(psb_out_unit,'("Solving matrix       : ",a)')mname
  write(psb_out_unit,'("Function to compute  : ",a)')fname
  write(psb_out_unit,'("Algorithmic variant  : ",a)')variant
  write(psb_out_unit,'("Scaling              : ",f8.2)')scaling

  ! Read matrix from file
  call mm_mat_read(a,info,iunit=iunit,filename=mname)
  n = a%get_nrows()
  nnz = a%get_nzeros()
  write(psb_out_unit,'("Matrix of size matrix ",i6,"x",i6," nnz =",i6)')n,n,nnz

  ! Set the options
  call fun%set("FNAME",fname,info)
  call fun%set("VARIANT",variant,info)
  call fun%set("SCALING",scaling,info)

  ! Allocate vectors for size
  allocate(x(n), stat=info)
  if (info /= 0) write(psb_err_unit,*) "x(n): Allocation request denied"
  allocate(y(n), stat=info)
  if (info /= 0) write(psb_err_unit,*) "y(n): Allocation request denied"

  ! Compute the matrix function on the all one vector
  x(1:n) = 1_psb_dpk_
  t1 = psb_wtime()
  call fun%apply(a,y,x,info)
  t2 = psb_wtime()

  write(psb_out_unit,'("Elapsed time: ",es12.5)')t2-t1

  ! Check if we have to write to file
  if(dump) then
    call mm_array_write(y, "serialtest out vector", info, iunit=iunit, filename="result.mtx")
  end if

  ! Clean the memory and close
  if (allocated(x)) deallocate(x, stat=info)
  if (info /= 0) write(psb_err_unit,*) "x(n): Deallocation request denied"
  if (allocated(y)) deallocate(y, stat=info)
  if (info /= 0) write(psb_err_unit,*) "y(n): Deallocation request denied"


  call psb_exit(ictxt)
  stop

9999 call psb_error(ictxt)

  stop

  contains

    subroutine  get_parms(ictxt,mname,fname,variant,scaling,dump)
      ! This subroutine reads the parameters needed to run the serialtest
      ! program from standard input
      integer(psb_ipk_), intent(in)  :: ictxt
      character(len=*),  intent(out) :: mname,fname,variant
      real(psb_dpk_), intent(out)    :: scaling
      logical, intent(out)           :: dump

      character(len=100)    :: dumpchar
      integer(psb_ipk_)   :: ip, inp_unit
      character(len=1024) :: filename

      if (command_argument_count()>0) then
        call get_command_argument(1,filename)
        inp_unit = 30
        open(inp_unit,file=filename,action='read',iostat=info)
        if (info /= 0) then
          write(psb_err_unit,*) 'Could not open file ',filename,' for input'
          call psb_abort(ictxt)
          stop
        else
          write(psb_err_unit,*) 'Opened file ',trim(filename),' for input'
        end if
      else
        inp_unit=psb_inp_unit
      end if
      read(inp_unit,*) ip
      if (ip == 5) then
        read(inp_unit,*) mname
        read(inp_unit,*) fname
        read(inp_unit,*) variant
        read(inp_unit,*) scaling
        read(inp_unit,*) dumpchar
        select case (psb_toupper(dumpchar))
          case ("T","TRUE")
            dump = .true.
          case default
            dump = .false.
        end select
      else
        ! wrong number of parameters on input: all hopes are lost
        call pr_usage(izero)
        call psb_abort(ictxt)
        stop 1
      end if

      if (inp_unit /= psb_inp_unit) then
        close(inp_unit)
      end if

      return

    end subroutine get_parms

    subroutine pr_usage(iout)
      ! Prints out information on incorrected program usage
      integer(psb_ipk_) :: iout
      write(iout,*)'incorrect parameter(s) found'
      write(iout,*)' usage:  serialtest mname fname variant &
           & scaling dump'
      write(iout,*)' where:'
      write(iout,*)'     mname:    MatrixMarket file name'
      write(iout,*)'     fname :   Matrix function name'
      write(iout,*)'     variant : Variant of the algorithm for the fname function'
      write(iout,*)'     scaling : Scalar scaling parameter  '
      write(iout,*)'     dump :    Save result on file T/F   '
    end subroutine pr_usage

end program serialtest
