
module simulationTime_class

  use numPrecision
  use genericProcedures, only : fatalError

  implicit none
  private

  !!
  !! Simple module to keep track of simulation time for time-dependent calculations
  !!
  !! Allows easy public access to time step and current time, and provides a single interface to
  !! change time step for calculations with a variable time step, rather than needing to update the
  !! time step separately in all required modules (source, material, etc.)
  !! 
  type, public :: simulationTime
    real(defReal) :: step = ONE
    real(defReal) :: now  = ZERO
  end type simulationTime

  type(simulationTime), public :: time

  public :: setStep
  public :: nextStep
  public :: timeStep
  public :: timeNow
  public :: timeLeft

contains

  !!
  !! Set time step
  !!
  subroutine setStep(dt)
    real(defReal), intent(in) :: dt
    character(100), parameter :: Here = 'set (timeStep_class.f90)'

    if (dt <= ZERO) call fatalError(Here, 'Time step must be positive')

    time % step = dt

  end subroutine setStep

  !!
  !! Advance time by one time step
  !!
  subroutine nextStep()

    time % now = time % now + time % step

  end subroutine nextStep

  !!
  !! Return time step size
  !!
  function timeStep() result(dt)
    real(defReal) :: dt

    dt = time % step

  end function timeStep

  function timeNow() result(t)
    real(defReal) :: t

    t = time % now

  end function timeNow

  !!
  !! Return time remaining until end of time step
  !!
  function timeLeft(t) result(remaining_t)
    real(defReal), intent(in) :: t
    real(defReal)             :: remaining_t

    remaining_t = time % now + time % step - t

  end function timeLeft

end module simulationTime_class
