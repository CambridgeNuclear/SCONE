module tallyClerk_inter

  use numPrecision
  use genericProcedures,     only : fatalError
  use particle_class,        only : particle, phaseCoord
  use particleDungeon_class, only : particleDungeon

  implicit none
  private

  ! List of codes for diffrent reports
  integer(shortInt),parameter,public :: inColl_CODE     = 1000 ,&
                                        outColl_CODE    = 1001 ,&
                                        path_CODE       = 1002 ,&
                                        trans_CODE      = 1003 ,&
                                        hist_CODE       = 1004 ,&
                                        cycleStart_CODE = 1005 ,&
                                        cycleEnd_CODE   = 1006


  !!
  !! Abstract interface for a single tallyClerk.
  !! It recives reports from the admin and processed them into scores and estimates.
  !! Its responsibilites are as follows:
  !! 1) Score some result by accepting a subset of all avalible reports
  !! 2) Display implementation determined measure of convergance (usually some variance)
  !! 3) Can return information about reports it requires
  !!
  type, public,abstract :: tallyClerk
    private

  contains
    !!
    procedure :: reportInColl
    procedure :: reportOutColl
    procedure :: reportPath
    procedure :: reportTrans
    procedure :: reportHist
    procedure :: reportCycleStart
    procedure :: reportCycleEnd

    procedure(validReports), deferred :: validReports
    procedure(display), deferred      :: display
    !procedure(print),deferred         :: print *** Interface for this procedure will be defined shortly

  end type tallyClerk

  abstract interface
    !!
    !! Returns array of codes that represent diffrent reports
    !!
    function validReports(self) result(validCodes)
      import :: tallyClerk ,&
                shortInt
      class(tallyClerk),intent(in)               :: self
      integer(shortInt),dimension(:),allocatable :: validCodes
    end function validReports

    !!
    !! Display convergance progress on the console
    !!
    subroutine display(self)
      import :: tallyClerk
      class(tallyClerk), intent(in)  :: self
    end subroutine display

  end interface

contains

  !!
  !! Process incoming collision report
  !!
  subroutine reportInColl(self,p)
    class(tallyClerk), intent(inout)      :: self
    class(particle), intent(in)           :: p
    character(100),parameter    :: Here = 'reportInColl (tallyClerk_inter.f90)'

    call fatalError(Here,'Report was send to an instance that does not support it.')

  end subroutine reportInColl


  !!
  !! Process outgoing collision report
  !!
  subroutine reportOutColl(self,pre,post,MT,muL)
    class(tallyClerk), intent(inout)      :: self
    class(phaseCoord), intent(in)         :: pre
    class(particle), intent(in)           :: post
    integer(shortInt), intent(in)         :: MT
    real(defReal), intent(in)             :: muL
    character(100),parameter  :: Here = 'reportOutColl (tallyClerk_inter.f90)'

    call fatalError(Here,'Report was send to an instance that does not support it.')

  end subroutine reportOutColl

  !!
  !! Process pathlength report
  !! ASSUMPTIONS:
  !! Pathlength must be contained within a single cell and material
  !!
  subroutine reportPath(self,pre,post,cellId,L)
    class(tallyClerk), intent(inout)     :: self
    class(phaseCoord), intent(in)        :: pre
    class(particle), intent(in)          :: post
    integer(shortInt), intent(in)        :: cellId
    real(defReal), intent(in)            :: L
    character(100),parameter  :: Here = 'reportPath (tallyClerk_inter.f90)'

    call fatalError(Here,'Report was send to an instance that does not support it.')

  end subroutine reportPath

  !!
  !! Process transition report
  !! ASSUMPTIONS:
  !! Transition must be a straight line
  !! Pre and Post direction is assumed the same (aligned with r_pre -> r_post vector)
  !!
  subroutine reportTrans(self,pre,post)
    class(tallyClerk), intent(inout)     :: self
    class(phaseCoord), intent(in)        :: pre
    class(particle), intent(in)          :: post
    character(100),parameter  :: Here = 'reportTrans (tallyClerk_inter.f90)'

    call fatalError(Here,'Report was send to an instance that does not support it.')

  end subroutine reportTrans

  !!
  !! Process history report
  !! ASSUMPTIONS:
  !! **** FATE CODES NEED TO BE SPECIFIED
  !!
  subroutine reportHist(self,pre,post,fate)
    class(tallyClerk), intent(inout) :: self
    class(phaseCoord), intent(in)        :: pre
    class(particle), intent(in)          :: post
    integer(shortInt),intent(in)         :: fate
    character(100),parameter  :: Here = 'reportHist (tallyClerk_inter.f90)'

    call fatalError(Here,'Report was send to an instance that does not support it.')

  end subroutine reportHist

  !!
  !! Process beggining of a cycle
  !!
  subroutine reportCycleStart(self,start)
    class(tallyClerk), intent(inout) :: self
    class(particleDungeon), intent(in)   :: start
    character(100),parameter  :: Here = 'reportCycleStart (tallyClerk_inter.f90)'

    call fatalError(Here,'Report was send to an instance that does not support it.')

  end subroutine reportCycleStart

  !!
  !! Process end of the cycle
  !!
  subroutine reportCycleEnd(self,end)
    class(tallyClerk), intent(inout) :: self
    class(particleDungeon), intent(in)   :: end
    character(100),parameter  :: Here = 'reportCycleEnd (tallyClerk_inter.f90)'

    call fatalError(Here,'Report was send to an instance that does not support it.')

  end subroutine reportCycleEnd

end module tallyClerk_inter
