module psfun_krylov_mod
    ! The psfun_krylov_mod contains all the operations relative to the
    ! construction of the various type of Krylov subspaces. The idea is that
    ! every type of subspace and orthogonalization algorithms corresponds to a
    ! new type. For each of this types the fundamental operations are the
    ! 1) creation/allocation : given a sparse matrix and a sparse vector of the
    !                          opportune type, and the maximum size of the
    !                          Krylov subspace the routine inizializes and
    !                          allocates the memory for the basis and the
    !                          projected matrix.
    ! 2) grow                : the method grow grows the space of the given
    !                          number of basis vectors applying the
    !                          orthogonalization algorithm, until the allocated
    !                          size is reached
    ! 3) free                : deallocates and frees the memory occupied by the
    !                          basis, and by the projection matrix. The method
    !                          does not touch nor the sparse matrix neither the
    !                          vector
    ! 4) apply               : Apply V(1:k) or transpose(V(1:k)) to a given
    !                          vector
    ! 5) getproj             : Gets the projection matrix : it is meant to be
    !                          used for computing f(v(:,k)^T a (v(:,k)) with
    !                          an opportune method

  use psb_base_mod

  type :: psfun_d_lanczos
    ! Type psfun_d_lanczos is the type encoding the Krylov subspace for the
    ! sparse marix a and the vectr v of maximum size maxsize and whose actual
    ! number of basis vector is k. At step k the following lanczos relation
    ! holds:
    !   v(1:k)' A v(1:k) =  t(1:k,1:k)
    ! with t(1:k,1:k) a symmetric tridiagonal matrix
    type(psb_dspmat_type)               :: a       ! Distributed sparse matrix
    type(psb_d_vect_type)               :: b       ! Distributed dense vector
    type(psb_d_vect_type), allocatable  :: v(:)    ! Krylov basis
    type(psb_desc_type)                 :: desc_a  ! PSBLAS communicator
    integer(psb_ipk_)                   :: maxsize ! Krylov maximum size
    integer(psb_ipk_)                   :: k       ! Krylov actual vector
    integer(psb_ipk_)                   :: reorth  ! Type of reorthogonalization
    real(psb_dpk_), allocatable         :: t(:,:)  ! Projected matrix
  contains
    ! Constructor:
    procedure, pass(ksp) :: new      => psfun_d_lanczos_new
    ! procedure, pass(ksp) :: grow     => psfun_d_lanczos_grow             todo
    ! Getters
    ! procedure, pass(ksp) :: getproj  => psfun_d_lanczos_getproj          todo
    ! Operations
    ! procedure, pass(ksp) :: apply    => psfun_d_lanczos_apply            todo
    ! Destroyer
    ! procedure, pass(ksp) :: free     => psfun_d_lanczos_free             todo
  end type psfun_d_lanczos

  ! The operations for the type are contained in the relative submodule that is
  ! contained in polynomial/d_lanczos.
  interface
    module subroutine psfun_d_lanczos_new(ksp,a,b,desc_a,info,maxsize,reorth)
      class(psfun_d_lanczos), intent(inout)   :: ksp
      type(psb_dspmat_type), intent(in)       :: a
      type(psb_d_vect_type), intent(in)       :: b
      type(psb_desc_type), intent(in)         :: desc_a
      integer(psb_ipk_), intent(out)          :: info
      integer(psb_ipk_), intent(in), optional :: maxsize
      integer(psb_ipk_), intent(in), optional :: reorth
    end subroutine
  end interface


end module psfun_krylov_mod
