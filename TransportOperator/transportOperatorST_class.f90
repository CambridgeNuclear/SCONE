!!
!! Transport operator for surface tracking
!!
module transportOperatorST_class
  use numPrecision
  use universalVariables

  use errors_mod,                 only : fatalError
  use particle_class,             only : particle
  use particleDungeon_class,      only : particleDungeon
  use dictionary_class,           only : dictionary

  ! Superclass
  use transportOperator_inter,    only : transportOperator, init_super => init

  ! Geometry interfaces
  use geometry_inter,             only : geometry, distCache

  ! Tally interface
  use tallyCodes
  use tallyAdmin_class,           only : tallyAdmin

  ! Nuclear data interfaces
  use nuclearDatabase_inter,      only : nuclearDatabase

  implicit none
  private

  !!
  !! Transport operator that moves a particle with surface tracking
  !!
  !! Sample Input Dictionary:
  !!   trans { type transportOperatorST; cache 0;}
  !!
  type, public, extends(transportOperator) :: transportOperatorST
    logical(defBool)  :: cache = .true.
  contains
    procedure :: transit => surfaceTracking
    ! Override procedure
    procedure :: init
  end type transportOperatorST

contains

  !!
  !! Performs surface tracking until a collision point is found
  !!
  subroutine surfaceTracking(self, p, tally, thisCycle, nextCycle)
    class(transportOperatorST), intent(inout) :: self
    class(particle), intent(inout)            :: p
    type(tallyAdmin), intent(inout)           :: tally
    class(particleDungeon),intent(inout)      :: thisCycle
    class(particleDungeon),intent(inout)      :: nextCycle
    integer(shortInt)                         :: event
    real(defReal)                             :: sigmaT, dist, sigmaTrack, invSigmaTrack
    type(distCache)                           :: cache
    real(defReal), parameter                  :: tol  = 1.0E-12
    character(100), parameter :: Here = 'surfaceTracking (transportOperatorST_class.f90)'

    STLoop: do
        
      ! Get local conditions
      p % T = self % geom % getTemperature(p % coords)
      p % rho  = self % geom % getDensity(p % coords)
      
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
        if (dist /= dist) then
          print *, "Particle location: ", p % rGlobal()
          print *, "Particle direction: ", p % dirGlobal()
          print *, "Total XS: ", sigmaT
          call fatalError(Here, "Distance is NaN")
        end if

      end if

      ! Save state before movement
      call p % savePrePath()

      ! Move to the next stop.
      if (self % cache) then
        call self % geom % move_withCache(p % coords, dist, event, cache)

      else
        call self % geom % move(p % coords, dist, event)

      end if

      ! Send tally report for a path moved
      call tally % reportPath(p, dist)

      select case(p % matIdx())
      
        ! Kill particle if it has leaked
        case(OUTSIDE_FILL)
          p % isDead = .true.
          p % fate = LEAK_FATE

        ! Give error if the particle somehow ended in an undefined material
        case(UNDEF_MAT)
          print *, "Particle location: ", p % rGlobal()
          call fatalError(Here, "Particle is in undefined material")
        
        ! Give error if the particle is in a region with overlapping cells
        case(OVERLAP_MAT)
          print *, "Particle location: ", p % rGlobal()
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
  !! Initialise ST operator from a dictionary
  !!
  !! See transportOperator_inter for details
  !!
  subroutine init(self, dict)
    class(transportOperatorST), intent(inout) :: self
    class(dictionary), intent(in)             :: dict

    ! Initialise superclass
    call init_super(self, dict)

    if (dict % isPresent('cache')) then
      call dict % get(self % cache, 'cache')
    end if

  end subroutine init

end module transportOperatorST_class
