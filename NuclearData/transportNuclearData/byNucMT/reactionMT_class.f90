module reactionMT_class

  use numPrecision
  use genericProcedures,  only : fatalError
  use RNG_class,          only : RNG
  use aceCard_class,      only : aceCard
  use emissionENDF_class, only : emissionENDF


  implicit none
  private

  !!
  !! Type that stores all data relevant to an MT reaction
  !!
  type, public :: reactionMT
    private
    integer(shortInt),public               :: MT
    integer(shortInt)                      :: firstIdx
    real(defReal),dimension(:),allocatable :: xs
    type(emissionENDF)                     :: kinematics

  contains
    ! Interface procedures
    procedure :: getXS
    procedure :: sampleAngleEnergy
    procedure :: isInCMframe
    procedure :: addOnGrid

    ! Build procedures
    procedure :: init => init_fromACE

  end type reactionMT

contains

  !!
  !! Returns xs given index in a grid and interpolation factor
  !!
  function getXS(self,idx,f) result(xs)
    class(reactionMT), intent(in) :: self
    integer(shortInt), intent(in) :: idx
    real(defReal),intent(in)      :: f
    real(defReal)                 :: xs
    integer(shortInt)             :: loc_idx

    ! Obtain local index
    loc_idx = idx - self % firstIdx + 1

    ! Return 0.0 if below first index. Interpolate xs otherwise
    if( loc_idx <= 0) then
      xs = ZERO

    else
      xs = (ONE-f) * self % xs(loc_idx) + f * self % xs(loc_idx+1)

    end if

  end function getXS

  !!
  !! Returns .true. if kinematics data for the reaction is in Centre-of-Mass frame
  !!
  function isInCMframe(self) result(isIt)
    class(reactionMT), intent(in) :: self
    logical(defBool)              :: isIt

    isIT = self % kinematics % isInCMframe()

  end function isInCMframe

  !!
  !! Sample angle and energy of outgoing particle given incident energy and random number generator
  !! Translates the call to embedded emissionENDF
  !!
  subroutine sampleAngleEnergy(self,mu,E_out,E_in,rand)
    class(reactionMT), intent(in) :: self
    real(defReal), intent(out)    :: mu
    real(defReal), intent(out)    :: E_out
    real(defReal), intent(in)     :: E_in
    class(RNG), intent(inout)     :: rand

    call self % kinematics % sampleAngleEnergy(mu,E_out,E_in,rand)

  end subroutine sampleAngleEnergy

  !!
  !! Adds internal XS on the supplied grid
  !!
  subroutine addOnGrid(self,grid)
    class(reactionMT), intent(in)            :: self
    real(defReal),dimension(:),intent(inout) :: grid
    character(100),parameter :: Here ='addOnGrid (reactionMT_class.f90)'

    ! Check if the provided grid is  too small
    if(size(grid) < size(self % xs) + self % firstIdx -1) then
      call fatalError(Here,'Trying to add xs on a grid that is too small to fit them')
    end if

    ! Add XS
    grid(self % firstIdx:) = grid(self % firstIdx:) + self % xs

  end subroutine addOnGrid

  !!
  !! Initialise reactionMT from ACE
  !! aceCard read head can be in any position
  !! Read head will be moved in subroutine
  !!
  subroutine init_fromACE(self,ACE,MT)
    class(reactionMT),intent(inout) :: self
    type(aceCard),intent(inout)     :: ACE
    integer(shortInt),intent(in)    :: MT
    integer(shortInt)               :: N

    ! Load MT number
    self % MT = MT

    ! Read First Index and Number of Points
    self % firstIdx = ACE % firstIdxMT(MT)
    N = ACE % numXsPointsMT(MT)

    ! Read XSs

    self % xs = ACE % xsMT(MT)

    ! Read Kinematic data
    call self % kinematics % init(ACE,MT)

  end subroutine init_fromACE

end module reactionMT_class
