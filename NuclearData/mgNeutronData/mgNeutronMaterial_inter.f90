module mgNeutronMaterial_inter

  use numPrecision
  use genericProcedures, only : fatalError
  use RNG_class,         only : RNG
  use particle_class,    only : particle

  ! Nuclear Data Handles
  use materialHandle_inter,    only : materialHandle
  use neutronMaterial_inter,   only : neutronMaterial
  use neutronXsPackages_class, only : neutronMacroXSs

  ! MG NEUTRON CACHE
  use mgNeutronCache_mod,      only : cache_materialCache => materialCache

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
  !! Private Members:
  !!   fissile -> flag set to .true. if material is fissile
  !!
  !! Interface:
  !!   materialHandle interface
  !!   neutroNMaterial interface
  !!   getMacroXSs -> Get macroscopic XSs directly from group number and RNG
  !!   set         -> Sets fissile flag
  !!
  type, public, abstract, extends(neutronMaterial) :: mgNeutronMaterial
    private
    logical(defBool) :: fissile = .false.

  contains
    ! Superclass procedures
    procedure :: kill
    generic   :: getMacroXSs => getMacroXSs_byG
    procedure :: getMacroXSs_byP

    ! Local procedures
    procedure(getMacroXSs_byG), deferred    :: getMacroXSs_byG
    procedure(getTotalXS), deferred         :: getTotalXS
    procedure                               :: isFissile
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
    subroutine getMacroXSs_byG(self, xss, G, rand)
      import :: mgNeutronMaterial, neutronMacroXSs, shortInt, RNG
      class(mgNeutronMaterial), intent(in) :: self
      type(neutronMacroXSs), intent(out)   :: xss
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
      import :: mgNeutronMaterial, defReal, shortInt, RNG
      class(mgNeutronMaterial), intent(in) :: self
      integer(shortInt), intent(in)        :: G
      class(RNG), intent(inout)            :: rand
      real(defReal)                        :: xs
    end function getTotalXS
  end interface



contains

  !!
  !! Return Macroscopic XSs for the material given particle
  !!
  !! See neutronMaterial_inter for details
  !!
  subroutine getMacroXSs_byP(self, xss, p)
    class(mgNeutronMaterial), intent(in) :: self
    type(neutronMacroXSs), intent(out)   :: xss
    class(particle), intent(in)          :: p
    character(100), parameter :: Here = 'getMacroXSs_byP (mgNeutronMateerial_inter.f90)'

    if (.not. p % isMG) call fatalError(Here, 'CE particle was given to MG data')

    associate (matCache => cache_materialCache(p % matIdx()))

      if (matCache % G_tail /= p % G) then
        ! Get cross sections
        call self % getMacroXSs(xss, p % G, p % pRNG)
        ! Update cache
        matCache % xss = xss
        matCache % G_tail = p % G

      else
        ! Retrieve cross sections from cache
        xss = matCache % xss

      end if

    end associate

  end subroutine getMacroXSs_byP

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
