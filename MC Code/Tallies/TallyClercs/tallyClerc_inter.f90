module tallyClerc_inter

  use numPrecision
  use genericProcedures,     only : fatalError
  use particle_class,        only : particle, phaseCoord
  use particleDungeon_class, only : particleDungeon

  implicit none
  private

  ! List of codes for diffrent reports
  integer(shortInt),parameter :: collision_CODE  = 1000 ,&
                                 path_CODE       = 1001 ,&
                                 trans_CODE      = 1002 ,&
                                 hist_CODE       = 1003 ,&
                                 cycleStart_CODE = 1004 ,&
                                 cycleEnd_CODE   = 1005


  !!
  !! Abstract interface for a single tallyClerc.
  !! It recives reports from the admin and processed them into scores and estimates.
  !! Its responsibilites are as follows:
  !! 1) Score some result by accepting a subset of all avalible reports
  !! 2) Display implementation determined measure of convergance (usually some variance)
  !! 3) Can return information about reports it requires
  !!
  type, public,abstract :: tallyClerc
    private

  contains
    !!
    procedure :: reportCollision
    procedure :: reportPath
    procedure :: reportTrans
    procedure :: reportHist
    procedure :: reportCycleStart
    procedure :: reportCycleEnd

    procedure(validReports), deferred :: validReports
    procedure(display), deferred      :: display
    !procedure(print),deferred         :: print *** Interface for this procedure will be defined shortly

  end type tallyClerc

  abstract interface
    !!
    !! Returns array of codes that represent diffrent reports
    !!
    function validReports(self) result(validCodes)
      import :: tallyClerc ,&
                shortInt
      class(tallyClerc),intent(in)               :: self
      integer(shortInt),dimension(:),allocatable :: validCodes
    end function validReports

    !!
    !! Display convergance progress on the console
    !!
    subroutine display(self)
      import :: tallyClerc
      class(tallyClerc), intent(in)  :: self
    end subroutine display

  end interface

contains

  !!
  !! Process collision report
  !!
  subroutine reportCollision(self,pre,post,MT,muL)
    class(tallyClerc), intent(inout)      :: self
    class(phaseCoord), intent(in)         :: pre
    class(particle), intent(in)           :: post
    integer(shortInt), intent(in)         :: MT
    real(defReal), intent(in)             :: muL
    character(100),parameter  :: Here = 'reportCollision (tallyClerc_inter.f90)'

    call fatalError(Here,'Report was send to an instance that does not support it.')

  end subroutine reportCollision

  !!
  !! Process pathlength report
  !! ASSUMPTIONS:
  !! Pathlength must be contained within a single cell and material
  !!
  subroutine reportPath(self,pre,post,cellId,L)
    class(tallyClerc), intent(inout)     :: self
    class(phaseCoord), intent(in)        :: pre
    class(particle), intent(in)          :: post
    integer(shortInt), intent(in)        :: cellId
    real(defReal), intent(in)            :: L
    character(100),parameter  :: Here = 'reportPath (tallyClerc_inter.f90)'

    call fatalError(Here,'Report was send to an instance that does not support it.')

  end subroutine reportPath

  !!
  !! Process transition report
  !! ASSUMPTIONS:
  !! Transition must be a straight line
  !! Pre and Post direction is assumed the same (aligned with r_pre -> r_post vector)
  !!
  subroutine reportTrans(self,pre,post)
    class(tallyClerc), intent(inout)     :: self
    class(phaseCoord), intent(in)        :: pre
    class(particle), intent(in)          :: post
    character(100),parameter  :: Here = 'reportTrans (tallyClerc_inter.f90)'

    call fatalError(Here,'Report was send to an instance that does not support it.')

  end subroutine reportTrans

  !!
  !! Process history report
  !! ASSUMPTIONS:
  !! **** FATE CODES NEED TO BE SPECIFIED
  !!
  subroutine reportHist(self,pre,post,fate)
    class(tallyClerc), intent(inout) :: self
    class(phaseCoord), intent(in)        :: pre
    class(particle), intent(in)          :: post
    integer(shortInt),intent(in)         :: fate
    character(100),parameter  :: Here = 'reportHist (tallyClerc_inter.f90)'

    call fatalError(Here,'Report was send to an instance that does not support it.')

  end subroutine reportHist

  !!
  !! Process beggining of a cycle
  !!
  subroutine reportCycleStart(self,start)
    class(tallyClerc), intent(inout) :: self
    class(particleDungeon), intent(in)   :: start
    character(100),parameter  :: Here = 'reportCycleStart (tallyClerc_inter.f90)'

    call fatalError(Here,'Report was send to an instance that does not support it.')

  end subroutine reportCycleStart

  !!
  !! Process end of the cycle
  !!
  subroutine reportCycleEnd(self,end)
    class(tallyClerc), intent(inout) :: self
    class(particleDungeon), intent(in)   :: end
    character(100),parameter  :: Here = 'reportCycleEnd (tallyClerc_inter.f90)'

    call fatalError(Here,'Report was send to an instance that does not support it.')

  end subroutine reportCycleEnd

end module tallyClerc_inter
