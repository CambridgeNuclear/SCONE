module tallyAdminBase_class

  use numPrecision
  use particle_class,        only : particle, phaseCoord
  use particleDungeon_class, only : particleDungeon

  implicit none
  private

  !!
  !! Base class for the tallies black box.
  !! Its responsibilities are as flolow:
  !! 1) Accept events reports and routes then to individual tallyClercs
  !! 2) Returns k-eff estimate for a current cycle
  !! 3) Controls end of calculation
  !! 4) Controls printing of calculation progress (to a console)
  !! 5) Controls printing of result estimators to a file (filePath and file  Format)
  !!
  !! This class will be extended by inheritance to provide additional functionality
  !! i.e. return mesh based weight windows based on fission matrix or similar
  !!
  !!
  type, public:: tallyAdminBase
  contains
    ! Report Interface
    procedure :: reportCollision
    procedure :: reportPath
    procedure :: reportTrans
    procedure :: reportHist
    procedure :: reportCycleStart
    procedure :: reportCycleEnd

    ! Access Procedures
    procedure :: k_eff

    ! Display procedures
    !procedure :: display

    ! File writing procedures
    !procedure :: print

    ! Build procedures
    !procedure :: init
  end type tallyAdminBase
    
contains

  !!
  !! Process collision report
  !!
  subroutine reportCollision(self,pre,post,MT,muL)
    class(tallyAdminBase), intent(inout)  :: self
    class(phaseCoord), intent(in)         :: pre
    class(particle), intent(in)           :: post
    integer(shortInt), intent(in)         :: MT
    real(defReal), intent(in)             :: muL
  end subroutine reportCollision

  !!
  !! Process pathlength report
  !! ASSUMPTIONS:
  !! Pathlength must be contained within a single cell and material
  !!
  subroutine reportPath(self,pre,post,cellId,L)
    class(tallyAdminBase), intent(inout) :: self
    class(phaseCoord), intent(in)        :: pre
    class(particle), intent(in)          :: post
    integer(shortInt), intent(in)        :: cellId
    real(defReal), intent(in)            :: L
  end subroutine reportPath

  !!
  !! Process transition report
  !! ASSUMPTIONS:
  !! Transition must be a straight line
  !! Pre and Post direction is assumed the same (aligned with r_pre -> r_post vector)
  !!
  subroutine reportTrans(self,pre,post)
    class(tallyAdminBase), intent(inout) :: self
    class(phaseCoord), intent(in)        :: pre
    class(particle), intent(in)          :: post
  end subroutine reportTrans

  !!
  !! Process history report
  !! ASSUMPTIONS:
  !! **** FATE CODES NEED TO BE SPECIFIED
  !!
  subroutine reportHist(self,pre,post,fate)
    class(tallyAdminBase), intent(inout) :: self
    class(phaseCoord), intent(in)        :: pre
    class(particle), intent(in)          :: post
    integer(shortInt),intent(in)         :: fate

  end subroutine reportHist

  !!
  !! Process beggining of a cycle
  !!
  subroutine reportCycleStart(self,start)
    class(tallyAdminBase), intent(inout) :: self
    class(particleDungeon), intent(in)   :: start

  end subroutine reportCycleStart

  !!
  !! Process end of the cycle
  !!
  subroutine reportCycleEnd(self,end)
    class(tallyAdminBase), intent(inout) :: self
    class(particleDungeon), intent(in)   :: end

  end subroutine reportCycleEnd


  !!
  !! Return estimate of k_eff
  !!
  subroutine k_eff(self,k,STD)
    class(tallyAdminBase), intent(in)   :: self
    real(defReal), intent(out)          :: k
    real(defReal),optional, intent(out) :: STD

  end subroutine k_eff

end module tallyAdminBase_class
