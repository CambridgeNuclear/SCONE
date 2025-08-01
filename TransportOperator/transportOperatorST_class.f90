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
    real(defReal)                             :: sigmaT, dist, speed, time
    type(distCache)                           :: cache
    character(100), parameter :: Here = 'surfaceTracking (transportOperatorST_class.f90)'

    
    STLoop: do

      ! Obtain the local cross-section
      if (p % matIdx() == VOID_MAT) then
        dist = INFINITY

      else
        sigmaT = self % xsData % getTrackingXS(p, p % matIdx(), MATERIAL_XS)
        dist = -log( p % pRNG % get()) / sigmaT

        ! Should never happen! Catches NaN distances
        if (dist /= dist) call fatalError(Here, "Distance is NaN")

        ! Set a max flight distance due to hitting the time-boundary
        if (p % timeMax > ZERO) then
          speed = p % getSpeed()
          time = dist / speed + p % time

          ! Hit the time boundary: set max distance and age the particle
          if (time > p % timeMax) then
            dist = speed * (p % timeMax - p % time)
            p % fate = AGED_FATE
          
          ! Advance particle in time
          else
            p % time = time
          end if

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

      ! Kill particle if it has leaked
      if (p % matIdx() == OUTSIDE_FILL) then
        p % isDead = .true.
        p % fate = LEAK_FATE
      end if

      ! Give error if the particle somehow ended in an undefined material
      if (p % matIdx() == UNDEF_MAT) then
        print *, p % rGlobal()
        call fatalError(Here, "Particle is in undefined material")
      end if

      ! Return if particle is stopped by collision, death, or aging
      if (event == COLL_EV .or. p % isDead .or. p % fate == AGED_FATE) exit STLoop

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
