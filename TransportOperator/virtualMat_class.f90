
module virtualMat_class

  use numPrecision
  use universalVariables

  use genericProcedures,          only : fatalError, numToChar

  use dynArray_class,             only : dynIntArray

  use nuclearDatabase_inter,      only : nuclearDatabase

  use particle_class,             only : particle, P_PHOTON, P_MATERIAL

  implicit None
  private

  !!
  !!
  !!
  type, public                      :: virtualMat
    type(dynIntArray)               :: realMats
    real(defReal)                   :: majorant_inv
    class(nuclearDatabase), pointer :: xsData
  contains
    procedure           :: init
    procedure           :: addRealMat
    procedure           :: updateMajorant
  end type virtualMat

contains

  subroutine init(self, nucData)
    class(virtualMat), intent(inout)            :: self
    class(nuclearDatabase), pointer, intent(in) :: nucData

    self % xsData => nucData

  end subroutine


  subroutine addRealMat(self, matIdx)
    class(virtualMat), intent(inout) :: self
    integer(shortInt), intent(in)    :: matIdx

    ! Add matIdx to array if not already present
    if (.not. self % realMats % isPresent(matIdx)) then
      call self % realMats % add(matIdx)
    end if

  end subroutine addRealMat


  subroutine updateMajorant(self)
    class(virtualMat), intent(inout) :: self
    integer(shortInt)                :: i
    real(defReal)                    :: majorant, sigma
    class(particle), allocatable     :: p
    character(100), parameter :: Here = ''

    allocate(p)
    !p % type = P_PHOTON
    !p % isMG = .true.
    ! TODO: Loop through groups when doing MG simulations
    p % G = 1

    majorant = ZERO

    ! Find majorant opacity of virtual material
    do i = 1, self % realMats % getSize()
      sigma = self % xsData % getTransMatXS(p, self % realMats % get(i))
      if (sigma > majorant) majorant = sigma
    end do

    self % majorant_inv = ONE/majorant

  end subroutine updateMajorant

end module virtualMat_class
