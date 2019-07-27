module datalessMaterials_class

  use numPrecision
  use genericProcedures, only : linFind, searchError, fatalError
  use dictionary_class,  only : dictionary
  use nuclearData_inter, only : nuclearData

  implicit none
  private


  !!
  !! Dummy nuclear data object
  !! Stores material names and indices without any nuclear data information
  !!
  type, public,extends(nuclearData) :: datalessMaterials
    private
    character(nameLen), dimension(:), allocatable :: materials ! the names of each material

  contains
    procedure :: init
    procedure :: kill
    procedure :: getIdx
    procedure :: getName

  end type datalessMaterials

contains

  !!
  !! Initialisation from a dictionary
  !!
  subroutine init(self, dict, matNames)
    class(datalessMaterials), intent(inout)     :: self
    class(dictionary), intent(in)               :: dict
    character(nameLen),dimension(:), intent(in) :: matNames

    if(allocated(self % materials)) deallocate(self % materials)

    self % materials = matNames

  end subroutine init

  !!
  !! Deallocate space
  !!
  elemental subroutine kill(self)
    class(datalessMaterials), intent(inout) :: self

    if(allocated(self % materials)) deallocate (self % materials)

  end subroutine kill


  !!
  !! Returns material index for given material name
  !! Throws error if material is not found
  !!
  function getIdx(self,matName) result (matIdx)
    class(datalessMaterials), intent(in)    :: self
    character(*), intent(in)                :: matName
    integer(shortInt)                       :: matIdx
    character(100), parameter               :: Here='getIdx (datalessMaterial_class.f90)'

    matIdx = linFind(self % materials, matName)
    call searchError(matIdx,Here)

  end function getIdx

  !!
  !! Returns material name for given material index
  !! Throws error if material index does not correspond to valid material
  !!
  function getName(self,matIdx) result(matName)
    class(datalessMaterials), intent(in)    :: self
    integer(shortInt), intent(in)           :: matIdx
    character(nameLen)                      :: matName
    character(100), parameter               :: Here='getName (datalessMaterials_class.f90)'

    if ( (matIdx < 1) .or. (matIdx > size(self % materials)) ) then
      call fatalError(Here,'Provided matIdx is invalid')

    end if

    matName = self % materials(matIdx)

  end function getName


end module datalessMaterials_class
