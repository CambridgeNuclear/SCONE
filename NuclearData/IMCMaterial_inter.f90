module IMCMaterial_inter

  use numPrecision
  use particle_class,       only : particle

  ! Nuclear Data Interfaces
  use materialHandle_inter,    only : materialHandle
  use IMCXsPackages_class, only : IMCMacroXSs

  implicit none
  private

  !!
  !! Public Pointer Cast
  !!
  public :: IMCMaterial_CptrCast

  !!
  !! Abstract interface far all IMC Materials (CE and MG)
  !!
  !! It was created to expose access to some key information in the context of
  !! tallying where one is not interested whether MG or CE data is used
  !!
  !! Interface:
  !!   materialHandle interface
  !!   isFissle    -> Return true if material is fissile
  !!   getMacroXSs -> Return Macroscopic XSs given particle with energy data
  !!
  type, public, abstract, extends(materialHandle) :: IMCMaterial
    private
  contains
    generic                              :: getMacroXSs => getMacroXSs_byP
    procedure(isFissile),       deferred :: isFissile
    procedure(getMacroXSs_byP), deferred :: getMacroXSs_byP
  end type IMCMaterial

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
      import :: IMCMaterial, defBool
      class(IMCMaterial), intent(in) :: self
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
      import :: IMCMaterial, particle, IMCMacroXSs
      class(IMCMaterial), intent(in) :: self
      type(IMCMacroXSs), intent(out) :: xss
      class(particle), intent(in)        :: p
    end subroutine getMacroXSs_byP

  end interface

contains


  !!
  !! Cast materialHandle pointer to IMCMaterial pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class materialHandle
  !!
  !! Result:
  !!   Null is source is not of IMCMaterial
  !!   Pointer to source if source is IMCMaterial class
  !!
  pure function IMCMaterial_CptrCast(source) result(ptr)
    class(materialHandle), pointer, intent(in) :: source
    class(IMCMaterial), pointer            :: ptr

    select type(source)
      class is(IMCMaterial)
        ptr => source

      class default
        ptr => null()
    end select

  end function IMCMaterial_CptrCast


end module IMCMaterial_inter