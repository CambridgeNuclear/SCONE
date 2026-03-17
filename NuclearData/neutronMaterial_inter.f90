module neutronMaterial_inter

  use numPrecision
  use particle_class,       only : particle

  ! Nuclear Data Interfaces
  use materialHandle_inter,    only : materialHandle
  use neutronXsPackages_class, only : neutronMacroXSs

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public :: neutronMaterial_CptrCast

  !!
  !! Abstract interface far all neutron Materials (CE and MG)
  !!
  !! It was created to expose access to some key information in the context of
  !! tallying where one is not interested whether MG or CE data is used
  !!
  !! Interface:
  !!   materialHandle interface
  !!   isFissle    -> Return true if material is fissile
  !!   getMacroXSs -> Return Macroscopic XSs given particle with energy data
  !!
  type, public, abstract, extends(materialHandle) :: neutronMaterial
    private
  contains
    generic                              :: getMacroXSs => getMacroXSs_byP
    procedure(isFissile),       deferred :: isFissile
    procedure(getMacroXSs_byP), deferred :: getMacroXSs_byP
    procedure(getMTxs), deferred         :: getMTxs

  end type neutronMaterial

  abstract interface
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
      import :: neutronMaterial, defBool
      class(neutronMaterial), intent(in) :: self
      logical(defBool)                   :: isIt
    end function isFissile

    !!
    !! Return Macroscopic XSs for the material given particle
    !!
    !! Args:
    !!   xss [out]    -> Cross section package to store the data
    !!   p [in]       -> Particle that provides energy or energy group
    !!
    !! Errors:
    !!   fatalError if energy value/group is outside bounds
    !!   fatalError if MG particle is given to CE data and vice versa
    !!
    subroutine getMacroXSs_byP(self, xss, p)
      import :: neutronMaterial, particle, neutronMacroXSs
      class(neutronMaterial), intent(in) :: self
      type(neutronMacroXSs), intent(out) :: xss
      class(particle), intent(in)        :: p
    end subroutine getMacroXSs_byP

    !!
    !! Return Macroscopic XS for the material given particle and an MT number
    !!
    !! Args:
    !!   MT [in]  -> MT number
    !!   p [in]   -> Particle that provides energy or energy group
    !!
    !! Errors:
    !!   fatalError if energy value is outside bounds
    !!   fatalError if MG particle is given to CE data
    !!
    !! NOTE: despite being in the interface, this function only makes sense
    !!       for CE. The MG extension returns a fatalError if called
    !!
    function getMTxs(self, MT, p) result(xs)
      import :: neutronMaterial, particle, shortInt, defReal
      class(neutronMaterial), intent(in) :: self
      integer(shortInt), intent(in)      :: MT
      class(particle), intent(in)        :: p
      real(defReal)                      :: xs
    end function getMTxs

  end interface

contains


  !!
  !! Cast materialHandle pointer to neutronMaterial pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class materialHandle
  !!
  !! Result:
  !!   Null is source is not of neutronMaterial
  !!   Pointer to source if source is neutronMaterial class
  !!
  pure function neutronMaterial_CptrCast(source) result(ptr)
    class(materialHandle), pointer, intent(in) :: source
    class(neutronMaterial), pointer            :: ptr

    select type(source)
      class is(neutronMaterial)
        ptr => source

      class default
        ptr => null()
    end select

  end function neutronMaterial_CptrCast


end module neutronMaterial_inter
