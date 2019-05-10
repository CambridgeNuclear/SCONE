module energyFilter_class

  use numPrecision
  use genericProcedures, only : fatalError, numToChar
  use particle_class,    only : particleState
  use dictionary_class,  only : dictionary
  use tallyFilter_inter, only : tallyFilter

  implicit none
  private

  !!
  !! Filter that tests if energy of a state is in a closed, energy interval.
  !! Returns .false. for MG particles
  !!
  !! NOTE: Energy is not enforced to be +ve. Could be useful in debugging
  !!
  type, public,extends(tallyFilter) :: energyFilter
    private
    real(defReal) :: Emin
    real(defReal) :: Emax
  contains
    procedure :: init
    procedure :: isPass

    !! Instance specific procedures
    procedure :: build

  end type energyFilter

contains

  !!
  !! Initialise energyFilter from dictionary
  !!
  subroutine init(self,dict)
    class(energyFilter), intent(inout) :: self
    class(dictionary), intent(in)      :: dict
    real(defReal)                      :: E1, E2

    ! Get bounds
    call dict % get(E1,'Emin')
    call dict % get(E2,'Emax')

    call self % build(E1, E2)

  end subroutine init

  !!
  !! Returns true if energy value is between specified bounds
  !!
  elemental function isPass(self,state) result(passed)
    class(energyFilter), intent(in)  :: self
    class(particleState), intent(in) :: state
    logical(defBool)                 :: passed
    real(defReal)                    :: E

    ! Special case - MG particle
    if(state % isMG) then
      passed = .false.
      return
    end if

    ! Normal case - CE paricle
    E = state % E
    passed = (self % Emin <= E) .and. (E <= self % Emax)

  end function isPass

  !!
  !! Build energyFilter from components
  !!
  !! Args:
  !!   Emin [in] -> minimum energy [MeV]
  !!   Emax [in] -> maximum energy [MeV]
  !!
  !! Errors:
  !!   fatalError if Emin > Emax
  !!
  subroutine build(self, Emin, Emax)
    class(energyFilter), intent(inout) :: self
    real(defReal), intent(in)          :: Emin
    real(defReal), intent(in)          :: Emax
    character(100), parameter :: Here = 'build (energyFilter_class.f90)'

    self % Emin = Emin
    self % Emax = Emax

    ! Verify bounds
    if( self % Emax <= self % Emin) then
      call fatalError(Here,'Emin='// numToChar(self % Emin) //' is larger or equal to Emax=' // numToChar(self % Emax))
    end if

  end subroutine build
    
end module energyFilter_class
