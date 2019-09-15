module mgNeutronMaterial_inter

  use numPrecision
  use RNG_class, only : RNG

  ! Nuclear Data Handles
  use materialHandle_inter,    only : materialHandle
  use neutronXsPackages_class, only : neutronMacroXSs

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public :: mgNeutronMaterial_CptrCast

  !!
  !! Extendable procedures is subclasses
  !!
  public :: kill

  !!
  !! Abstract interface for all MG neutron Materials
  !!
  !!
  !! Interface:
  !!   materialHandle interface
  !!   isFissile -> return .true. if is a fissile material
  !!   set       -> Sets fissile flag
  !!
  type, public, abstract, extends(materialHandle) :: mgNeutronMaterial
    logical(defBool) :: fissile = .false.

  contains
    ! Superclass procedures
    procedure :: kill

    ! Local procedures
    procedure(getMacroXSs), deferred        :: getMacroXSs
    procedure, non_overridable              :: isFissile
    procedure                               :: set

  end type mgNeutronMaterial



  abstract interface

    !!
    !! Return Macroscopic XSs for the material
    !!
    !! Args:
    !!   xss [out]    -> Cross section package to store the data
    !!   G [in]       -> Requested energy group
    !!   rand [inout] -> Random Number Generator
    !!
    !! Errors:
    !!   fatalError if G is out-of-bounds for the stored data
    !!
    subroutine getMacroXSs(self, xss, G, rand)
      import :: mgNeutronMaterial, neutronMacroXSs, shortInt, RNG
      class(mgNeutronMaterial), intent(in) :: self
      type(neutronMacroXSs), intent(out)   :: xss
      integer(shortInt), intent(in)        :: G
      class(RNG), intent(inout)            :: rand
    end subroutine getMacroXSs


  end interface



contains

  !!
  !! Return .true. if the MG material is fissile
  !!
  !! Args:
  !!   None
  !!
  !! Errors:
  !!   None
  !!
  elemental function isFissile(self) result(isIt)
    class(mgNeutronMaterial), intent(in) :: self
    logical(defBool)                     :: isIt

    isIt = self % fissile

  end function isFissile

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(mgNeutronMaterial), intent(inout) :: self

    self % fissile = .false.

  end subroutine kill

  !!
  !! Set fissile flag
  !!
  !! All arguments are optional. Use with keyword association e.g.
  !!   call mat % set( fissile = .true.)
  !!
  !! Args:
  !!   fissile [in] -> flag indicating whether fission data is present
  !!
  subroutine set(self, fissile)
    class(mgNeutronMaterial), intent(inout) :: self
    logical(defBool), intent(in), optional  :: fissile

    if(present(fissile)) self % fissile = fissile

  end subroutine set

  !!
  !! Cast materialHandle pointer to mgNeutronMaterial pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class materialHandle
  !!
  !! Result:
  !!   Null is source is not of ceNeutronMaterial
  !!   Pointer to source if source is ceNeutronMaterial class
  !!
  pure function mgNeutronMaterial_CptrCast(source) result(ptr)
    class(materialHandle), pointer, intent(in) :: source
    class(mgNeutronMaterial), pointer          :: ptr

    select type(source)
      class is(mgNeutronMaterial)
        ptr => source

      class default
        ptr => null()
    end select

  end function mgNeutronMaterial_CptrCast

end module mgNeutronMaterial_inter
