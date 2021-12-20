module mgIMCMaterial_inter

  use numPrecision
  use genericProcedures, only : fatalError
  use RNG_class,         only : RNG
  use particle_class,    only : particle

  ! Nuclear Data Handles
  use materialHandle_inter,    only : materialHandle
  use IMCMaterial_inter,       only : IMCMaterial
  use IMCXsPackages_class,     only : IMCMacroXSs

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public :: mgIMCMaterial_CptrCast

  !!
  !! Extendable procedures is subclasses
  !!
  public :: kill

  !!
  !! Abstract interface for all MG IMC Materials
  !!
  !! Interface:
  !!   materialHandle interface
  !!   neutroNMaterial interface
  !!   getMacroXSs -> Get macroscopic XSs directly from group number and RNG
  !!
  type, public, abstract, extends(IMCMaterial) :: mgIMCMaterial
    private

  contains
    ! Superclass procedures
    procedure :: kill
    generic   :: getMacroXSs => getMacroXSs_byG
    procedure :: getMacroXSs_byP

    ! Local procedures
    procedure(getMacroXSs_byG), deferred    :: getMacroXSs_byG
    procedure(getTotalXS), deferred         :: getTotalXS

  end type mgIMCMaterial



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
    subroutine getMacroXSs_byG(self, xss, G, rand)
      import :: mgIMCMaterial, IMCMacroXSs, shortInt, RNG
      class(mgIMCMaterial), intent(in)     :: self
      type(IMCMacroXSs), intent(out)       :: xss
      integer(shortInt), intent(in)        :: G
      class(RNG), intent(inout)            :: rand
    end subroutine getMacroXSs_byG

    !!
    !! Return Macroscopic Total XS for the material
    !!
    !! Args:
    !!   G [in]       -> Requested energygroup
    !!   rand [inout] -> Random number generator
    !!
    !! Errors:
    !!   fatalError if G is out-of-bounds for the stored data
    !!
    function getTotalXS(self, G, rand) result(xs)
      import :: mgIMCMaterial, defReal, shortInt, RNG
      class(mgIMCMaterial), intent(in)     :: self
      integer(shortInt), intent(in)        :: G
      class(RNG), intent(inout)            :: rand
      real(defReal)                        :: xs
    end function getTotalXS
  end interface



contains

  !!
  !! Return Macroscopic XSs for the material given particle
  !!
  !! See IMCMaterial_inter for details
  !!
  subroutine getMacroXSs_byP(self, xss, p)
    class(mgIMCMaterial), intent(in)     :: self
    type(IMCMacroXSs), intent(out)       :: xss
    class(particle), intent(in)          :: p
    character(100), parameter :: Here = 'getMacroXSs_byP (mgIMCMateerial_inter.f90)'

    if( p % isMG) then
      call self % getMacroXSs(xss, p % G, p % pRNG)

    else
      call fatalError(Here, 'CE particle was given to MG data')

    end if
  end subroutine getMacroXSs_byP

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(mgIMCMaterial), intent(inout) :: self

  end subroutine kill

  !!
  !! Cast materialHandle pointer to mgIMCMaterial pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class materialHandle
  !!
  !! Result:
  !!   Null is source is not of ceIMCMaterial
  !!   Pointer to source if source is ceIMCMaterial class
  !!
  pure function mgIMCMaterial_CptrCast(source) result(ptr)
    class(materialHandle), pointer, intent(in) :: source
    class(mgIMCMaterial), pointer              :: ptr

    select type(source)
      class is(mgIMCMaterial)
        ptr => source

      class default
        ptr => null()
    end select

  end function mgIMCMaterial_CptrCast

end module mgIMCMaterial_inter
