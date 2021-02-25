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
!
! This file is based on the getp file from:
!
!                Parallel Sparse BLAS  version 3.5
!      (C) Copyright 2006-2018
!        Salvatore Filippone
!        Alfredo Buttari
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions
!    are met:
!      1. Redistributions of source code must retain the above copyright
!         notice, this list of conditions and the following disclaimer.
!      2. Redistributions in binary form must reproduce the above copyright
!         notice, this list of conditions, and the following disclaimer in the
!         documentation and/or other materials provided with the distribution.
!      3. The name of the PSBLAS group or the names of its contributors may
!         not be used to endorse or promote products derived from this
!         software without specific written permission.
!
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!    ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!    PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
!    BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!    POSSIBILITY OF SUCH DAMAGE.
!
!
Module getp
  interface get_parms
    module procedure get_dparms
  end interface

contains
  !
  ! Get iteration parameters from the command line
  !
  subroutine  get_dparms(ctxt,mtrx_file,rhs_file,eta,zeta,filefmt,kmethd,ptype,part,&
       & afmt,istopc,itmax,itrace,irst,eps)
    use psb_base_mod
    type(psb_ctxt_type) :: ctxt
    character(len=2)  :: filefmt
    character(len=40) :: kmethd, mtrx_file, rhs_file, ptype
    character(len=20) :: part
    integer(psb_ipk_) :: iret, istopc,itmax,itrace,irst
    character(len=40) :: charbuf
    real(psb_dpk_) :: eta,zeta,eps
    character    :: afmt*5
    integer(psb_ipk_) :: np, iam
    integer(psb_ipk_) :: inparms(40), ip, inp_unit
    character(len=1024)   :: filename

    call psb_info(ctxt,iam,np)
    if (iam == 0) then
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
      ! Read Input Parameters
      read(inp_unit,*) ip
      if (ip >= 9) then
        read(inp_unit,*) mtrx_file
        read(inp_unit,*) rhs_file
        read(inp_unit,*) eta
        read(inp_unit,*) zeta
        read(inp_unit,*) filefmt
        read(inp_unit,*) kmethd
        read(inp_unit,*) ptype
        read(inp_unit,*) afmt
        read(inp_unit,*) part


        call psb_bcast(ctxt,mtrx_file)
        call psb_bcast(ctxt,rhs_file)
        call psb_bcast(ctxt,eta)
        call psb_bcast(ctxt,zeta)
        call psb_bcast(ctxt,filefmt)
        call psb_bcast(ctxt,kmethd)
        call psb_bcast(ctxt,ptype)
        call psb_bcast(ctxt,afmt)
        call psb_bcast(ctxt,part)

        if (ip >= 10) then
          read(inp_unit,*) istopc
        else
          istopc=1
        endif
        if (ip >= 11) then
          read(inp_unit,*) itmax
        else
          itmax=500
        endif
        if (ip >= 12) then
          read(inp_unit,*) itrace
        else
          itrace=-1
        endif
        if (ip >= 13) then
          read(inp_unit,*) irst
        else
          irst  = 1
        endif
        if (ip >= 14) then
          read(inp_unit,*) eps
        else
          eps=1.d-6
        endif
        inparms(1) = istopc
        inparms(2) = itmax
        inparms(3) = itrace
        inparms(4) = irst
        call psb_bcast(ctxt,inparms(1:4))
        call psb_bcast(ctxt,eps)

        write(psb_out_unit,'("Solving matrix       : ",a)')  mtrx_file
        write(psb_out_unit,'("Number of processors : ",i3)') np
        write(psb_out_unit,'("Data distribution    : ",a)') part
        write(psb_out_unit,'("Iterative method     : ",a)')  kmethd
        write(psb_out_unit,'("Preconditioner       : ",a)')  ptype
        write(psb_out_unit,'("Restart parameter    : ",i2)') irst
        write(psb_out_unit,'("Storage format       : ",a)')  afmt(1:3)
        write(psb_out_unit,'(" ")')
      else
        write(psb_err_unit,*) 'Wrong format for input file'
        call psb_abort(ctxt)
        stop 1
      end if
      if (inp_unit /= psb_inp_unit) then
        close(inp_unit)
      end if
    else
      ! Receive Parameters
      call psb_bcast(ctxt,mtrx_file)
      call psb_bcast(ctxt,rhs_file)
      call psb_bcast(ctxt,eta)
      call psb_bcast(ctxt,zeta)
      call psb_bcast(ctxt,filefmt)
      call psb_bcast(ctxt,kmethd)
      call psb_bcast(ctxt,ptype)
      call psb_bcast(ctxt,afmt)
      call psb_bcast(ctxt,part)

      call psb_bcast(ctxt,inparms(1:4))
      istopc =  inparms(1)
      itmax  =  inparms(2)
      itrace =  inparms(3)
      irst   =  inparms(4)
      call psb_bcast(ctxt,eps)

    end if

  end subroutine get_dparms

end module getp
