module physicsPackage_inter

  use numPrecision
  use dictionary_class, only : dictionary

  implicit none
  private

  !!
  !! Abstract interface of physics Package
  !! Physics package is controles a calculation flow
  !! Each type of calculation has different physics package
  !! Loud is for displaying calculation progress
  !!
  type, public,abstract :: physicsPackage
    private
    logical(defBool), public :: loud
  contains
    procedure(init), deferred :: init
    procedure(run),deferred   :: run
  end type physicsPackage

  abstract interface

    !!
    !! Initialise Physics Package from dictionary
    !!
    subroutine init(self,dict,loud)
      import :: physicsPackage, &
                dictionary, &
                defBool
      class(physicsPackage), intent(inout)   :: self
      class(dictionary), intent(inout)       :: dict
      logical(defBool), intent(in), optional :: loud
    end subroutine init

    !!
    !! Run calculation in the physics package
    !!
    subroutine run(self)
      import :: physicsPackage
      class(physicsPackage), intent(inout) :: self
    end subroutine run

    !!
    !! Deallocate memory used by physicsPackage
    !!
    subroutine kill(self)
      import :: physicsPackage
      class(physicsPackage), intent(inout) :: self
    end subroutine kill

  end interface
contains

end module physicsPackage_inter
