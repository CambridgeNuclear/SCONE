module nuclearData_inter

  use numPrecision
  use dictionary_class, only : dictionary

  implicit none
  private


  !!
  !! Base interface for nuclear Data Module
  !!
  !! It is neccesery to allow other modules to be oblivious to whether nuclear Data is MG or CE
  !! or it conforms to some other, more specific interface.
  !! but still be able to access stuff common to all nuclear data instances. Namly:
  !! -> Material Names
  !! -> Material Indixes
  !!
  type,abstract, public :: nuclearData
    private
  contains
    procedure(init), deferred     :: init
    procedure(getIdx), deferred   :: getIdx
    procedure(getName), deferred  :: getName
    procedure(kill), deferred     :: kill

  end type nuclearData

  abstract interface

    !!
    !! Initialise nuclear data
    !! from material dictionary
    !! and an array of material names, that specifies material indices
    !!
    subroutine init(self, dict, matNames)
      import :: nuclearData, &
                dictionary, &
                nameLen
      class(nuclearData), intent(inout)            :: self
      class(dictionary), intent(in)                :: dict
      character(nameLen), dimension(:), intent(in) :: matNames
    end subroutine init

    !!
    !! Returns material index for given material name
    !! Throws error if material is not found
    !!
    function getIdx(self,matName) result(matIdx)
      import :: nuclearData, &
                shortInt
      class(nuclearData), intent(in) :: self
      character(*),intent(in)        :: matName
      integer(shortInt)              :: matIdx

    end function getIdx

    !!
    !! Returns material name for given material index
    !! Throws error if material index does not correspond to valid material
    !!
    function getName(self,matIdx) result(matName)
      import :: nuclearData, &
                shortInt, &
                nameLen
      class(nuclearData), intent(in) :: self
      integer(shortInt), intent(in)  :: matIdx
      character(nameLen)             :: matName

    end function getName

    !!
    !! Deallocate memory occupied by nuclearData
    !!
    elemental subroutine kill(self)
      import :: nuclearData
      class(nuclearData), intent(inout) :: self

    end subroutine kill


  end interface
    
end module nuclearData_inter
