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
    procedure :: filter

  end type energyFilter

contains

  !!
  !! Initialise energyFilter from dictionary
  !!
  subroutine init(self,dict)
    class(energyFilter), intent(inout) :: self
    class(dictionary), intent(in)      :: dict
    real(defReal)                      :: temp
    character(100), parameter :: Here = 'init (energyFilter_class.f90)'

    ! Get bounds
    call dict % get(temp,'Emin')
    self % Emin = temp

    call dict % get(temp,'Emax')
    self % Emax = temp

    ! Verify bounds
    if( self % Emax <= self % Emin) then
      call fatalError(Here,'Emin='// numToChar(self % Emin) //' is larger or equal to Emax=' // numToChar(self % Emax))
    end if

  end subroutine init

  !!
  !! Returns true if energy value is between specified bounds
  !!
  elemental function filter(self,state) result(passed)
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

  end function filter
    
end module energyFilter_class
