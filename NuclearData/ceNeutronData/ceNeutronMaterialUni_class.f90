module ceNeutronMaterialUni_class

  use numPrecision
  use genericProcedures, only : fatalError, binarySearch, numToChar

  ! Nuclear Data Handles
  use ceNeutronMaterial_class,  only : ceNeutronMaterial
  use materialHandle_inter,     only : materialHandle

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public ceNeutronMaterialUni_TptrCast

  !!
  !! Extension of ceNeutronMaterial, used to write and read material-wise unionised energy grids.
  !! For more information see ceNeutronMaterial_class
  !!
  !! Called only by aceNeutronDatabaseUni and aceNeutronDatabaseUniIdx
  !!
  !! Procedures:
  !!   search -> return index and interpolation factor for the energy of interest in
  !!             the unionised energy grid
  !!
  type, public, extends(ceNeutronMaterial) :: ceNeutronMaterialUni
    real(defReal), dimension(:), allocatable     :: unionGrid
    real(defReal), dimension(:), allocatable     :: totalXS
    integer(shortInt),dimension(:,:),allocatable :: gridIdx
  contains
    procedure,non_overridable :: search
  end type ceNeutronMaterialUni

contains

  !!
  !! Search energy for grid and interpolation factor for energy E
  !!
  !! Interpolation factor definition:
  !!   f = (E - E_low) / (E_top - E_low)
  !!   E = E_top * f + E_low * (1-f)
  !!
  !! Args:
  !!   idx [out] -> index of the bottom bin for energy E
  !!   f   [out] -> value of the interpolation factor for energy E
  !!   E   [in]  -> Energy to search for [MeV]
  !!
  !! Fatal Errors:
  !!   If energy E is beyond range terminate with fatalError
  !!
  subroutine search(self, idx, f, E)
    class(ceNeutronMaterialUni), intent(in) :: self
    integer(shortInt), intent(out)          :: idx
    real(defReal), intent(out)              :: f
    real(defReal), intent(in)               :: E
    character(100), parameter :: Here = 'search (ceNeutronMaterialUni_class.f90)'

    idx = binarySearch(self % unionGrid, E)
    if(idx <= 0) then
      call fatalError(Here,'Failed to find energy: '//numToChar(E)//&
                           ' for material '//numToChar(self % matIdx))
    end if

    associate(E_top => self % unionGrid(idx + 1), E_low  => self % unionGrid(idx))
      f = (E - E_low) / (E_top - E_low)
    end associate

  end subroutine search

  !!
  !! Cast materialHandle pointer to ceNeutronMaterialUni pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class materialHandle
  !!
  !! Result:
  !!   Null is source is not of ceNeutronMaterialUni
  !!   Pointer to source if source is ceNeutronMaterialUni class
  !!
  pure function ceNeutronMaterialUni_TptrCast(source) result(ptr)
    class(materialHandle), pointer, intent(in) :: source
    type(ceNeutronMaterialUni), pointer        :: ptr

    select type(source)
    type is(ceNeutronMaterialUni)
        ptr => source

      class default
        ptr => null()
    end select

  end function ceNeutronMaterialUni_TptrCast


end module ceNeutronMaterialUni_class
