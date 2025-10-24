!!
!! Transport operator for hybrid tracking
!!
module transportOperatorHT_class
  use numPrecision
  use universalVariables

  use errors_mod,                 only : fatalError
  use genericProcedures,          only : numToChar
  use particle_class,             only : particle
  use particleDungeon_class,      only : particleDungeon
  use dictionary_class,           only : dictionary

  ! Tally interface
  use tallyCodes
  use tallyAdmin_class,           only : tallyAdmin

  ! Superclass
  use transportOperator_inter,    only : transportOperator, init_super => init

  ! Geometry interfaces
  use geometry_inter,             only : geometry

  ! Nuclear data interfaces
  use nuclearDataReg_mod,         only : ndReg_get => get
  use nuclearDatabase_inter,      only : nuclearDatabase

  implicit none
  private

  !!
  !! Transport operator that moves a particle with hybrid tracking
  !!
  type, public, extends(transportOperator) :: transportOperatorHT
    real(defReal) :: cutoff   ! Cutoff threshold between ST and DT

  contains
    procedure :: transit => tracking_selection
    procedure, private :: deltaTracking
    procedure, private :: surfaceTracking
    ! Override procedure
    procedure :: init

  end type transportOperatorHT

contains

  subroutine tracking_selection(self, p, tally, thisCycle, nextCycle)
    class(transportOperatorHT), intent(inout)              :: self
    class(particle), intent(inout)                         :: p
    type(tallyAdmin), intent(inout)                        :: tally
    class(particleDungeon), intent(inout)                  :: thisCycle
    class(particleDungeon), intent(inout)                  :: nextCycle
    real(defReal)                                          :: majorant_inv, sigmaT, ratio
    character(100), parameter :: Here = 'hybridTracking (transportOperatorHT_class.f90)'

    ! Get majornat XS inverse: 1/Sigma_majorant
    majorant_inv = ONE / self % xsData % getTrackingXS(p, p % matIdx(), MAJORANT_XS)

    ! Obtain the local cross-section. Always choose ST in void
    if (p % matIdx() == VOID_MAT) then
      sigmaT = ZERO
    else
      sigmaT = self % xsData % getTrackMatXS(p, p % matIdx())
    end if

    ! Calculate ratio between local cross-section and majorant
    ratio = sigmaT * majorant_inv

    ! Cut-off criterion to decide on tracking method
    if (ratio > (ONE - self % cutoff)) then
      call deltaTracking(self, p, tally, thisCycle, nextCycle)
    else
      call surfaceTracking(self, p, tally, thisCycle, nextCycle)
    end if

  end subroutine tracking_selection

  !!
  !! Performs delta tracking until a real collision point is found
  !!
  subroutine deltaTracking(self, p, tally, thisCycle, nextCycle)
    class(transportOperatorHT), intent(inout) :: self
    class(particle), intent(inout)            :: p
    type(tallyAdmin), intent(inout)           :: tally
    class(particleDungeon), intent(inout)     :: thisCycle
    class(particleDungeon), intent(inout)     :: nextCycle
    real(defReal)                             :: majorant_inv, sigmaT, distance
    character(100), parameter :: Here = 'deltaTracking (transportOperatorHT_class.f90)'

    ! Get majorant XS inverse: 1/Sigma_majorant
    majorant_inv = ONE / self % xsData % getTrackingXS(p, p % matIdx(), MAJORANT_XS)

   ! Should never happen! Prevents Inf distances
    if (abs(majorant_inv) > huge(majorant_inv)) call fatalError(Here, "Majorant is 0")

    DTLoop:do
      distance = -log( p % pRNG % get() ) * majorant_inv

      ! Move particle in the geometry
      call self % geom % teleport(p % coords, distance)

      select case(p % matIdx())

        ! If particle has leaked exit
        case(OUTSIDE_FILL)
          p % fate = LEAK_FATE
          p % isDead = .true.
          return


        ! Check for void
        case(VOID_MAT)
          call tally % reportInColl(p, .true.)
          cycle DTLoop

        ! Give error if the particle somehow ended in an undefined material
        case(UNDEF_MAT)
          print *, "Particle location: ", p % rGlobal()
          call fatalError(Here, "Particle is in undefined material")

        ! Give error if the particle somehow ended in an overlap material
        case(OVERLAP_MAT)
          print *, "Particle location: ", p % rGlobal()
          call fatalError(Here, "Particle is in overlapping cells")

        case default
          ! All is well

      end select

      ! Obtain the local cross-section
      sigmaT = self % xsData % getTrackMatXS(p, p % matIdx())

      ! Roll RNG to determine if the collision is real or virtual
      ! Exit the loop if the collision is real, report collision if virtual
      if (p % pRNG % get() < sigmaT*majorant_inv) then
        exit DTLoop
      else
        call tally % reportInColl(p, .true.)
      end if

    end do DTLoop

    call tally % reportTrans(p)

  end subroutine deltaTracking

  !!
  !! Performs surface tracking until a collision point is found
  !!
  subroutine surfaceTracking(self, p, tally, thisCycle, nextCycle)
    class(transportOperatorHT), intent(inout) :: self
    class(particle), intent(inout)            :: p
    type(tallyAdmin), intent(inout)           :: tally
    class(particleDungeon),intent(inout)      :: thisCycle
    class(particleDungeon),intent(inout)      :: nextCycle
    integer(shortInt)                         :: event
    real(defReal)                             :: sigmaT, dist, sigmaTrack, invSigmaTrack
    real(defReal), parameter                  :: tol  = 1.0E-12
    character(100), parameter :: Here = 'surfaceTracking (transportOperatorHT_class.f90)'

    STLoop: do
      
      sigmaTrack = self % xsData % getTrackingXS(p, p % matIdx(), MATERIAL_XS)

      ! Obtain the local cross-section, depending on the material
      ! This branch is called in the case of voids with no imposed XS
      if (sigmaTrack < tol) then
        
        dist = INFINITY
        invSigmaTrack = INFINITY
        sigmaT = ZERO

      else
      
        invSigmaTrack = ONE / sigmaTrack
        dist = -log( p % pRNG % get()) * invSigmaTrack
        
        ! Obtain the local cross-section
        sigmaT = self % xsData % getTrackMatXS(p, p % matIdx())

        ! Should never happen! Catches NaN distances
        if (dist /= dist) call fatalError(Here, "Distance is NaN")

      end if
      
      ! Save state before movement
      call p % savePrePath()

      ! Move to the next stop. NOTE: "move" resets dist to distanced moved!
      call self % geom % move(p % coords, dist, event)

      ! Send tally report for a path moved
      call tally % reportPath(p, dist)
      
      select case(p % matIdx())
      
        ! Kill particle if it has leaked
        case(OUTSIDE_FILL)
          p % isDead = .true.
          p % fate = LEAK_FATE

        ! Give error if the particle somehow ended in an undefined material
        case(UNDEF_MAT)
          print*, 'Particle location: ', p % rGlobal()
          call fatalError(Here, "Particle is in undefined material")

        ! Give error if the particle ended in an overlap material
        case(OVERLAP_MAT)
          print*, 'Particle location: ', p % rGlobal()
          call fatalError(Here, "Particle is in overlapping cells")
      
        case default
          ! All is well

      end select
      
      if (p % isDead) exit STLoop

      ! Roll RNG to determine if the collision is real or virtual
      ! Exit the loop if the collision is real, report collision if virtual
      if (event == COLL_EV) then
        if (p % pRNG % get() < sigmaT*invSigmaTrack) then
          exit STLoop
        else
          call tally % reportInColl(p, .true.)
        end if
      end if

    end do STLoop

    call tally % reportTrans(p)

  end subroutine surfaceTracking

  !!
  !! Initialise HT operator from a dictionary
  !!
  !! See transportOperator_inter for more details
  !!
  subroutine init(self, dict)
    class(transportOperatorHT), intent(inout) :: self
    class(dictionary), intent(in)             :: dict

    ! Initialise superclass
    call init_super(self, dict)

    ! Retrieve DT-ST probability cutoff
    call dict % getOrDefault(self % cutoff,'cutoff',0.9_defReal)

  end subroutine init


end module transportOperatorHT_class
