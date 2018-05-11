module tally_inter

  use numPrecision
  use particle_class, only : particle, phaseCoord

  implicit none
  private

  type, public,abstract :: tally
  contains
    ! Calculation contron (?) !*reconsider name
    procedure(startCycle),deferred      :: startCycle
    procedure(endCycle),deferred        :: endCycle

    ! Print Current Results
    procedure(printEstimates),deferred  :: printEstimates

    ! Report Interface
    procedure(reportCollision),deferred :: reportCollision
    procedure(reportPath),deferred      :: reportPath
    procedure(reportTrans),deferred     :: reportTrans
    procedure(reportHist),deferred      :: reportHist
  end type tally
    
  abstract interface

    !!
    !! Inform tally about a start of new calcultion cycle (generation)
    !! Provide starting particle population
    !!
    subroutine startCycle(self,startPop)
      import :: tally, &
                shortInt
      class(tally), intent(inout)   :: self
      integer(shortInt), intent(in) :: startPop
    end subroutine startCycle

    !!
    !! Inform tally about the end of calculation cycle (generation)
    !! Provide end particle population
    !! Obtain current best estimate of criticality.
    !! Best estimate DOES NOT have to be analog. It is specified by an implementation
    !! DOES NOT reset startPop (multiple endCycles for one startPop are allowed) *** MAY CHANGE
    !!
    subroutine endCycle(self,k_est,endPop)
      import :: tally, &
                defReal ,&
                shortInt
      class(tally), intent(inout)   :: self
      real(defReal),intent(out)     :: k_est
      integer(shortInt), intent(in) :: endPop
    end subroutine endCycle

    !!
    !! Prints scores to the console
    !! What scores are printed is configuaration dependant
    !!
    subroutine printScores(self)
      import :: tally
      class(tally), intent(inout)  :: self
    end subroutine printScores

    !!
    !! Send collision report
    !!
    subroutine reportCollision(self,pre,post,MT,muL)
      import :: tally, &
                particle ,&
                phaseCoord ,&
                defReal, &
                shortInt
      class(tally), intent(inout)   :: self
      class(phaseCoord), intent(in) :: pre
      class(particle), intent(in)   :: post
      integer(shortInt), intent(in) :: MT
      real(defReal), intent(in)     :: muL
    end subroutine reportCollision

    !!
    !! Send pathlengths report
    !! Pathlength must be contained within a single cell and material
    !!
    subroutine reportPath(self,pre,post,cellId,L)
      import :: tally, &
                particle ,&
                phaseCoord ,&
                defReal, &
                shortInt
      class(tally), intent(inout)   :: self
      class(phaseCoord), intent(in) :: pre
      class(particle), intent(in)   :: post
      integer(shortInt), intent(in) :: cellId
      real(defReal), intent(in)     :: L
    end subroutine reportPath

    !!
    !! Send transition report
    !! Transition must be a straight line
    !! Pre and Post direction is assumed the same (aligned with r_pre -> r_post vector)
    !!
    subroutine reportTrans(self,pre,post)
      import :: tally, &
                particle ,&
                phaseCoord ,&
                defReal
      class(tally), intent(inout)   :: self
      class(phaseCoord), intent(in) :: pre
      class(particle), intent(in)   :: post
    end subroutine reportTrans

    !!
    !! Send history report
    !!
    subroutine reportHist(self,pre,post)
      import :: tally, &
                particle ,&
                phaseCoord ,&
                defReal, &
                shortInt
      class(tally), intent(inout)   :: self
      class(phaseCoord), intent(in) :: pre
      class(particle), intent(in)   :: post
    end subroutine reportHist

    subroutine printEstimates(self)
      import :: tally
      class(tally), intent(inout)  :: self
    end subroutine printEstimates

  end interface

end module tally_inter
