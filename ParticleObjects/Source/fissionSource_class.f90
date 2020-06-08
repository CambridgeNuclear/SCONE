module fissionSource_class

  use numPrecision
  use universalVariables,      only: OUTSIDE_MAT, NOT_FOUND
  use genericProcedures,       only: fatalError
  use dictionary_class,        only: dictionary
  use RNG_class,               only: RNG

  use particle_class,          only: particleState, P_NEUTRON, P_PHOTON
  use source_inter,            only: source

  use geometry_inter,          only: geometry
  use materialMenu_mod,        only: matIdx
  use nuclearDatabase_inter
  use fissionCE_class,         only: fissionCE, fissionCE_TptrCast
  use fissionMG_class,         only: fissionMG, fissionMG_TptrCast
  use ceNeutronDatabase_inter, only: ceNeutronDatabase, ceNeutronDatabase_CptrCast
  use reactionHandle_inter,    only: reactionHandle
  use uncorrelatedReactionCE_inter, only: uncorrelatedReactionCE, uncorrelatedReactionCE_CptrCast
  use reactionMG_inter,        only: reactionMG_CptrCast, reactionMG

  implicit none
  private

  !!
  !! Neutron Source from distributed fission sites
  !!
  !! Private members:
  !!   isMG         -> is the source multi-group?
  !!
  !! Interface:
  !!   source_inter Interface
  !!
  type, public,extends(source) :: fissionSource
    private
    logical(defBool)            :: isMG   = .false.
    real(defReal), dimension(3) :: bottom = ZERO
    real(defReal), dimension(3) :: top = ZERO
  contains
    procedure :: init
    procedure :: sampleType
    procedure :: samplePosition
    procedure :: sampleEnergy
    procedure :: sampleEnergyAngle
    procedure :: kill
  end type fissionSource

contains

  !!
  !! Initialise general source
  !!
  !! Read dictionary to obtain source information and provide
  !! geometry to allow basic check or position sampling
  !!
  !! Args:
  !!   dict [in] -> dict containing point source information
  !!   geom [in] -> pointer to geometry, for checking that point source is inside geometry
  !!                or for sampling position if mat/cell are specified
  !!
  !! Result:
  !!   An initialised general source
  !!
  !! Errors:
  !!   - error if an unrecognised particle type is provided
  !!   - error if either direction or position have more than 3 components
  !!   - error if both CE and MG is specified
  !!   - error if position/cell/mat are specified with each other
  !!   - error if neither energy type is specified
  !!   - error if energy and reaction are specified
  !!
  subroutine init(self, dict, geom)
    class(fissionSource), intent(inout)      :: self
    class(dictionary), intent(in)            :: dict
    class(geometry), pointer, intent(in)     :: geom
    character(30)                            :: type
    character(100), parameter :: Here = 'init (fissionSource_class.f90)'

    ! ! Provide geometry info to source
    ! self % geom => geom
    !
    ! ! Must inform of energy type - assume CE otherwise
    ! call dict % getOrDefault(dataType,'data','ce')
    ! select case(dataType)
    !
    !   case('ce')
    !     self % isMG = .false.
    !
    !   case('mg')
    !     self % isMG = .true.
    !
    !   case default
    !     call fatalError(Here, 'Invalid source data type specified: must be ce or mg')
    ! end select

  end subroutine init

  !!
  !! Sample particle's phase space co-ordinates
  !!
  !! See source_inter for details
  !!
  subroutine sampleParticle(self, p, rand)
    class(fissionSource), intent(inout)  :: self
    type(particleState), intent(inout)   :: p
    class(RNG), intent(inout)            :: rand


  end subroutine sampleParticle

  !!
  !! Return to uninitialised state
  !!
  subroutine kill(self)
    class(fissionSource), intent(inout) :: self
    !
    ! self % geom => null()
    ! if allocated(self % bb) deallocate(self % bb)
    ! if allocated(self % r) deallocate(self % r)
    ! if allocated(self % dir) deallocate(self % dir)
    ! if associated(self % reacCE) self % reacCE => null()
    ! if associated(self % reacMG) self % reacMG => null()

  end subroutine kill

end module fissionSource_class
