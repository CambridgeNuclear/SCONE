module P1MG_class

  use numPrecision
  use endfConstants
  use genericProcedures,  only : fatalError, numToChar
  use legendrePoly_func,  only : sampleLegendre
  use RNG_class,          only : RNG
  use dictionary_class,   only : dictionary
  use IOdictionary_class, only : IOdictionary
  use isotropicMG_class,  only : isotropicMG, init_super         => init,        &
                                              readMaterial_super => readMaterial,&
                                              activeIdx_super    => activeIdx
  implicit none
  private

  !!
  !! Extension of isotropicMG class, which uses P1 scattering matrix when sampling
  !! the outgoing scattering angle
  !!
  type, public,extends(isotropicMG) :: P1MG
    private
    real(defReal),dimension(:,:,:),allocatable :: P1  ! (G_out, G_in, matIdx)
  contains
    ! Overwrite key isotropicMG procedures
    procedure :: init
    procedure :: readMaterial
    procedure :: sampleMuGout
  end type P1MG

contains

  !!
  !! Initialise P1MG
  !!
  subroutine init(self, dict)
    class(P1MG), intent(inout)    :: self
    class(dictionary), intent(in) :: dict
    integer(shortInt)             :: nG, nMat
    character(nameLen),dimension(:),allocatable :: matNames

    ! Read material names
    call dict % keysDict(matNames)

    ! Read number of energy groups and materials
    call dict % get(nG, 'numberOfGroups')
    nMat = size( matNames )

    ! Allocate space for normalised P1 coefficients
    allocate(self % P1(nG, nG, nMat))

    ! Call superclass initialisation (note that it will call extended readMaterial in this module)
    call init_super(self,dict)

  end subroutine init

  !!
  !! Extend superclass subroutine to read additional data to produce P1 transport corrected
  !! transport XSs
  !!
  subroutine readMaterial(self, dict, idx, nG)
    class(P1MG), intent(inout)           :: self
    class(dictionary), intent(in)          :: dict
    integer(shortInt), intent(in)          :: idx
    integer(shortInt), intent(in)          :: nG
    type(IOdictionary)                     :: xsDict
    integer(shortInt)                      :: i
    real(defReal),dimension(:),allocatable :: P0
    real(defReal),dimension(:),allocatable :: P1
    character(100), parameter :: Here='readMaterial (P1MG_class.f90)'

    ! Call superclass procedure
    call readMaterial_super(self, dict, idx, nG)

    ! Obtain path from dict and read into xsDict
    call xsDict % initFrom( dict % getChar('xsFile'))

    ! Obtain P0 and P1 scattering XSs
    call xsDict % get(P1, 'scattering_P1')
    call xsDict % get(P0, 'scattering_P0')

    ! Verify size of the P0 and P1 matrix
    if (size(P1) /= nG*nG) then
      call fatalError(Here,'scattering_P1 matrix contains only:' // &
                            numToChar(size(P1))// ' elements')
    else if (size(P0) /= nG*nG) then
      call fatalError(Here,'scattering_P0 matrix contains only:' // &
                            numToChar(size(P0))// ' elements')
    end if

    ! Obtain normalised Legandre P1 coefficiant
    do i=1,size(P0)
      if(P0(i) == ZERO) then
        P1(i) = ZERO

      else
        P1(i) = P1(i) / P0(i)

      end if
    end do

    ! Verify that values do not go -ve

    ! Reshape and store
    self % P1(:, :, idx) = reshape(P1, [nG, nG])

  end subroutine readMaterial

  !!
  !! Samples deflection angle in LAB frame and post-colission energy group
  !! This implementation is CRAP. SHOULD BE BRANCHLESS. IMPROVE IT! **##~~##**
  !!
  subroutine sampleMuGout(self,mu,G_out,G_in,rand,MT,matIdx)
    class(P1MG), intent(in)        :: self
    real(defReal), intent(out)     :: mu
    integer(shortInt), intent(out) :: G_out
    integer(shortInt), intent(in)  :: G_in
    class(RNG), intent(inout)      :: rand
    integer(shortInt), intent(in)  :: MT
    integer(shortInt), intent(in)  :: matIdx
    real(defReal)                  :: r1, r2
    character(100), parameter      :: Here ='sampleMuGout (P1MG_class.f90)'

    r1 = rand % get()

    select case(MT)
      case(macroAllScatter)
        G_out = self % P_transferMatrix(G_in, matIdx) % invert(r1)

        ! Sample outgoing mu value
        mu = sampleLegendre(self % P1(G_out, G_in, matIdx), rand)

      case(macroFission)
        G_out = self % P_chiValues(matIdx) % invert(r1)
        mu = TWO * rand % get() - ONE

      case default
        call fatalError(Here,'Unrecoginsed MT number')

    end select

  end subroutine sampleMuGout

    
end module P1MG_class
