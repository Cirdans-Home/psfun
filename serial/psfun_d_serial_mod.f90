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
module psfun_d_serial_mod
    ! This module contains the generic interfaces for the computation of the
    ! different matrix functions included in the library. The idea is that this
    ! modules computes, in a serial way, :math:`y = f(\alpha A)x`.

    use psb_base_mod
    implicit none

    type, public :: psfun_d_serial
        character(len=20)   :: fname   = 'EXP'
        character(len=20)   :: variant = 'EXPOKIT'
        real(psb_dpk_)      :: scaling = 1.0_psb_dpk_
        integer(psb_ipk_)   :: padedegree = 6_psb_ipk_
    contains
        ! Set the options
        procedure, pass(fun) :: setstring  => psfun_d_setstring
        procedure, pass(fun) :: setreal    => psfun_d_setreal
        procedure, pass(fun) :: setinteger => psfun_d_setinteger
        generic, public :: set => setstring, setreal, setinteger
        ! Computes the function
        procedure, pass(fun) :: apply      => psfun_d_serial_apply
    end type psfun_d_serial

    private :: psfun_d_setstring, psfun_d_setreal, psfun_d_setinteger

contains

    subroutine psfun_d_setstring(fun,what,val,info)
        ! Set function for setting options defined by a string
        use psb_base_mod
        implicit none

        class(psfun_d_serial), intent(inout) :: fun   ! Function object
        character(len=*), intent(in)         :: what  ! String of option to set
        character(len=*), intent(in)         :: val   ! Value of the string
        integer(psb_ipk_), intent(out)       :: info  ! Output flag

        info = psb_success_
        select case (psb_toupper(what))
        case ("FNAME")
            fun%fname = val
        case ("VARIANT")
            fun%variant = val
        case default
            info = psb_err_invalid_args_combination_
        end select

    end subroutine psfun_d_setstring

    subroutine psfun_d_setreal(fun,what,val,info)
        ! Set function for setting options defined by a real
        use psb_base_mod
        implicit none

        class(psfun_d_serial), intent(inout) :: fun   ! Function object
        character(len=*), intent(in)         :: what  ! String of option to set
        real(psb_dpk_), intent(in)           :: val   ! Real Value of the option
        integer(psb_ipk_), intent(out)       :: info  ! Output flag

        info = psb_success_
        select case (psb_toupper(what))
        case("SCALING")
            fun%scaling = val
        case default
            info = psb_err_invalid_args_combination_
        end select

    end subroutine psfun_d_setreal

    subroutine psfun_d_setinteger(fun,what,val,info)
        ! Set function for setting options defined by an integer
        use psb_base_mod
        implicit none

        class(psfun_d_serial), intent(inout) :: fun   ! Function object
        character(len=*), intent(in)         :: what  ! String of option to set
        integer(psb_dpk_), intent(in)        :: val   ! Integer Value of the option
        integer(psb_ipk_), intent(out)       :: info  ! Output flag

        info = psb_success_
        select case (psb_toupper(what))
        case("PADE_DEGREE")
            fun%padedegree = val
        case default
            info = psb_err_invalid_args_combination_
        end select

    end subroutine psfun_d_setinteger

    subroutine psfun_d_serial_apply(fun,a,y,x,info)
        ! This is the core of the function apply on a serial matrix to compute
        ! :math:`y = f(\alpha*A) x`. It calls on the specific routines
        ! implementing the different functions. It is the function to modify if
        ! ones want to interface a new function that was not previously
        ! available or a new algorithm (variant) for an already existing
        ! function.
        use psb_base_mod
        implicit none

        class(psfun_d_serial), intent(inout) :: fun ! Function information
        real(psb_dpk_), intent(in)           :: a(:,:) ! Matrix
        real(psb_dpk_), intent(in)           :: x(:) ! Input vector
        real(psb_dpk_), intent(out)          :: y(:) ! Output vector
        integer(psb_ipk_), intent(out)       :: info ! Information on the output

        ! local variables
        integer(psb_ipk_)               :: n,m,lwsp,iexph,ns,shapes(2)
        real(psb_dpk_), allocatable     :: fA(:,:),wsp(:)
        integer(psb_ipk_), allocatable  :: ipiv(:)
        ! local constants
        integer(psb_ipk_)       :: err_act
        character(len=20)       :: name

        info = psb_success_
        name = 'psb_krylov'
        call psb_erractionsave(err_act)


        n = size(a,1)
        m = size(a,2)
        if( n /= m) then
            info = psb_err_from_subroutine_
            call psb_errpush(info,name,a_err=trim(fun%fname))
            goto 9999
        end if
        shapes(1) = n
        shapes(2) = m

        select case (fun%fname)
        case ('EXP')
            select case (fun%variant)
            case ('TAYLOR')
                allocate(fA(n,m), stat=info)
                if ( info /= 0) then
                    call psb_errpush(info,name,a_err=trim(fun%variant))
                    goto 9999
                end if

                call r8mat_expm2( n, fun%scaling*a, fA )
                y = matmul(fA,x)

                if (allocated(fA)) deallocate(fA, stat=info)
                if ( info /= 0) then
                    call psb_errpush(info,name,a_err=trim(fun%variant))
                    goto 9999
                end if
            case ('MATLAB')
                allocate(fA(n,m), stat=info)
                if ( info /= 0) then
                    call psb_errpush(info,name,a_err=trim(fun%variant))
                    goto 9999
                end if

                call r8mat_expm1( n, fun%scaling*a, fA )
                y = matmul(fA,x)

                if (allocated(fA)) deallocate(fA, stat=info)
                if ( info /= 0) then
                    call psb_errpush(info,name,a_err=trim(fun%variant))
                    goto 9999
                end if
            case('GENPADE')
                lwsp = 4*n*m+fun%padedegree+1
                allocate(wsp(lwsp), stat=info)
                if (info /= 0) then
                    call psb_errpush(info,name,a_err=trim(fun%variant))
                    goto 9999
                end if
                allocate(ipiv(m), stat=info)
                if (info /= 0) then
                    call psb_errpush(info,name,a_err=trim(fun%variant))
                    goto 9999
                end if

                call DGPADM(fun%padedegree,n,fun%scaling,a,m,wsp,lwsp,ipiv,iexph,ns,info)
                if (info /= 0) then
                    call psb_errpush(info,name,a_err=trim(fun%variant))
                    goto 9999
                end if

                allocate(fA(n,m), stat=info)
                if ( info /= 0) then
                    call psb_errpush(info,name,a_err=trim(fun%variant))
                    goto 9999
                end if

                fA = reshape(wsp, shape=shapes)
                y = matmul(fA,x)

                if (allocated(wsp)) deallocate(wsp, stat=info)
                if (info /= 0) then
                    call psb_errpush(info,name,a_err=trim(fun%variant))
                    goto 9999
                end if
                if (allocated(ipiv)) deallocate(ipiv, stat=info)
                if (info /= 0) then
                    call psb_errpush(info,name,a_err=trim(fun%variant))
                    goto 9999
                end if
                if (allocated(fA)) deallocate(fA, stat=info)
                if ( info /= 0) then
                    call psb_errpush(info,name,a_err=trim(fun%variant))
                    goto 9999
                end if
            case default
                info = psb_err_from_subroutine_
                call psb_errpush(info,name,a_err=trim(fun%variant))
                goto 9999
            end select
        case default
            info = psb_err_from_subroutine_
            call psb_errpush(info,name,a_err=trim(fun%fname))
            goto 9999
        end select

        if(info /= psb_success_) then
          info = psb_err_from_subroutine_
          call psb_errpush(info,name,a_err=trim(fun%fname))
          goto 9999
        end if

        call psb_erractionrestore(err_act)
        return

        9999 continue

        return

    end subroutine

end module psfun_d_serial_mod
