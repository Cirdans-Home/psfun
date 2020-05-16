submodule (psfun_krylov_mod) d_lanczos

contains
  module subroutine psfun_d_lanczos_new(ksp,a,b,desc_a,info,maxsize,reorth)
    use psb_base_mod
    implicit none

    class(psfun_d_lanczos), intent(inout)   :: ksp
    type(psb_dspmat_type), intent(in)       :: a
    type(psb_d_vect_type), intent(in)       :: b
    type(psb_desc_type), intent(in)         :: desc_a
    integer(psb_ipk_), intent(out)          :: info
    integer(psb_ipk_), intent(in), optional :: maxsize
    integer(psb_ipk_), intent(in), optional :: reorth
    ! local variables
    real(psb_dpk_), allocatable         :: t(:,:)
    type(psb_d_vect_type), allocatable  :: v(:)

    info = -1

    allocate(t(maxsize+1,maxsize+1), stat=info)
    if ( info /= 0) then
      write(*,*) "t(maxsize,maxsize): Allocation request denied"
      ! Error handling
    end if

    call psb_geall(v,desc_a,info,n=maxsize+1)

    if ( info /= 0) then
      write(*,*) "v(maxsize): Allocation request denied"
      ! Error handling
    end if

    ! Assign input variables to the object
    ksp%a = a
    ksp%b = b
    ksp%t = t
    if(present(reorth)) then
      if ((reorth > 1).or.(reorth < 0)) then
        ksp%reorth = 0
      else
        ksp%reorth = reorth
      end if
    else
      ksp%reorth = 0
    endif
    if(present(maxsize)) then
      if(maxsize < 0) then
        ksp%maxsize = min(100,a%get_ncols())
      else
        ksp%maxsize = maxsize
      end if
    else
      ksp%maxsize = min(100,a%get_ncols())
    endif
    ksp%v = v
    ksp%k = 0

  end subroutine
end submodule d_lanczos
