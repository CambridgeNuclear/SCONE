module virtualMat_class

  use numPrecision
  use universalVariables

  use genericProcedures,          only : fatalError, numToChar

  use dynArray_class,             only : dynIntArray

  use nuclearDatabase_inter,      only : nuclearDatabase

  use particle_class,             only : particle, P_PHOTON, P_MATERIAL

  use materialEquations,          only : mgEnergyGrid

  implicit None
  private

  !!
  !! Fake material class used for multiGeom delta tracking. Fills upper geometry cells, and stores
  !! information about the real materials overlapped in the lower geometry, specifically which
  !! materials these are and the majorant_inv within the upper geometry cell.
  !!
  !! See transportOperatorGeomHT_class.f90 for example of use.
  !!
  type, public :: virtualMat
    type(dynIntArray)                        :: realMats
    real(defReal), dimension(:), allocatable :: majorant_inv
    class(nuclearDatabase), pointer          :: xsData
  contains
    procedure           :: init
    procedure           :: addRealMat
    procedure           :: updateMajorant
  end type virtualMat

contains

  !!
  !! Allocate space for majorants
  !!
  subroutine init(self, nucData)
    class(virtualMat), intent(inout)            :: self
    class(nuclearDatabase), pointer, intent(in) :: nucData
    integer(shortInt)                           :: nG

    self % xsData => nucData

    ! Allocate space for majorants
    ! TODO: Currently this method of obtaining nG is specific to MG IMC
    if(associated(mgEnergyGrid)) then
      nG = mgEnergyGrid % getSize()
    else
      nG = 1
    end if

    allocate(self % majorant_inv(nG))

  end subroutine

  !!
  !! Add matIdx to the list of real materials in the underlying geometry that overlap with the
  !! region occupied by this virtual material.
  !!
  subroutine addRealMat(self, matIdx)
    class(virtualMat), intent(inout) :: self
    integer(shortInt), intent(in)    :: matIdx

    ! Add matIdx to array if not already present
    if (.not. self % realMats % isPresent(matIdx)) then
      call self % realMats % add(matIdx)
    end if

  end subroutine addRealMat

  !!
  !! Update the majorant_inv for each group
  !!
  subroutine updateMajorant(self)
    class(virtualMat), intent(inout) :: self
    integer(shortInt)                :: i, G
    real(defReal)                    :: majorant, sigma
    class(particle), allocatable     :: p
    character(100), parameter :: Here = ''

    ! Create particle for call to obtain XS - no data needed except p % G
    allocate(p)
    p % isMG = .true.

    majorant = ZERO

    ! Find majorant opacity of virtual material for each group
    do G = 1, size(self % majorant_inv)
      p % G = G
      do i = 1, self % realMats % getSize()
        sigma = self % xsData % getTransMatXS(p, self % realMats % get(i))
        if (sigma > majorant) majorant = sigma
      end do

      self % majorant_inv(G) = ONE/majorant

    end do

  end subroutine updateMajorant

end module virtualMat_class
