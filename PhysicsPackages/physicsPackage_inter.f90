module physicsPackage_inter

  use numPrecision
  use dictionary_class, only : dictionary

  implicit none
  private

  !!
  !! Abstract interface of physics Package
  !! Physics package is controles a calculation flow
  !! Each type of calculation has diffrent physics package
  !!
  type, public,abstract :: physicsPackage
    private
  contains
    procedure(init), deferred :: init
    procedure(run),deferred   :: run
  end type physicsPackage

  abstract interface

    !!
    !! Initialise Physics Package from dictionary
    !!
    subroutine init(self,dict)
      import :: physicsPackage, &
                dictionary
      class(physicsPackage), intent(inout) :: self
      class(dictionary), intent(inout)     :: dict
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
