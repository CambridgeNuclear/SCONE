module IMCMaterial_inter

  use numPrecision
  use particle_class,       only : particle
  use RNG_class,            only : RNG

  ! Nuclear Data Interfaces
  use materialHandle_inter,    only : materialHandle
  use IMCXsPackages_class,     only : IMCMacroXSs

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
  !!   getMacroXSs -> Return Macroscopic XSs given particle with energy data
  !!
  type, public, abstract, extends(materialHandle) :: IMCMaterial
    private
  contains
    generic                              :: getMacroXSs => getMacroXSs_byP
    procedure(getMacroXSs_byP), deferred :: getMacroXSs_byP
    procedure(updateMat), deferred       :: updateMat
    procedure(getEmittedRad), deferred   :: getEmittedRad
    procedure(getFleck), deferred        :: getFleck
    procedure(getTemp), deferred         :: getTemp
    procedure(setTimeStep), deferred     :: setTimeStep
  end type IMCMaterial

  abstract interface

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
      class(IMCMaterial), intent(in)     :: self
      type(IMCMacroXSs), intent(out)     :: xss
      class(particle), intent(in)        :: p
    end subroutine getMacroXSs_byP

    !!
    !! Update material properties at each time step
    !! First update energy using simple balance, then solve for temperature,
    !!  then update temperature-dependent properties
    !!
    !! Args:
    !!   tallyEnergy [in] -> Energy absorbed into material
    !!   printUpdate [in, optional] -> Bool, if true then will print updates to screen
    !!
    subroutine updateMat(self, tallyEnergy, printUpdate)
      import :: IMCMaterial, defReal, defBool
      class(IMCMaterial), intent(inout)      :: self
      real(defReal), intent(in)              :: tallyEnergy
      logical(defBool), intent(in), optional :: printUpdate
    end subroutine updateMat

    !!
    !! Return the equilibrium radiation energy density, U_r
    !!
    function getEmittedRad(self) result(emittedRad)
      import :: IMCMaterial, defReal, RNG
      class(IMCMaterial), intent(inout)  :: self
      real(defReal)                      :: emittedRad
    end function getEmittedRad

    !!
    !! Get Fleck factor of material
    !!
    function getFleck(self) result(fleck)
      import :: IMCMaterial, defReal
      class(IMCMaterial), intent(in) :: self
      real(defReal)                  :: fleck
    end function getFleck

    !!
    !! Get temperature of material
    !!
    function getTemp(self) result(T)
      import :: IMCMaterial, defReal
      class(IMCMaterial), intent(inout) :: self
      real(defReal)                     :: T
    end function getTemp

    !!
    !! Provide material with time step size
    !!
    !! Args:
    !!   dt [in] -> time step size [s]
    !!
    subroutine setTimeStep(self, dt)
      import :: IMCMaterial, defReal
      class(IMCMaterial), intent(inout) :: self
      real(defReal), intent(in)         :: dt
    end subroutine setTimeStep

  end interface

contains


  !!
  !! Cast materialHandle pointer to IMCMaterial pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class materialHandle
  !!
  !! Result:
  !!   Null if source is not of IMCMaterial
  !!   Pointer to source if source is IMCMaterial class
  !!
  pure function IMCMaterial_CptrCast(source) result(ptr)
    class(materialHandle), pointer, intent(in) :: source
    class(IMCMaterial), pointer                :: ptr

    select type(source)
      class is(IMCMaterial)
        ptr => source

      class default
        ptr => null()
    end select

  end function IMCMaterial_CptrCast


end module IMCMaterial_inter
