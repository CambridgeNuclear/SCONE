!!
!! Transport operator for delta tracking
!!
module transportOperatorDT_class
  use numPrecision
  use universalVariables

  use errors_mod,                 only : fatalError
  use genericProcedures,          only : numToChar
  use particle_class,             only : particle
  use particleDungeon_class,      only : particleDungeon
  use dictionary_class,           only : dictionary

  ! Superclass
  use transportOperator_inter,    only : transportOperator, init_super => init

  ! Geometry interfaces
  use geometry_inter,             only : geometry

  ! Tally interface
  use tallyCodes
  use tallyAdmin_class,           only : tallyAdmin

  ! Nuclear data interfaces
  use nuclearDataReg_mod,         only : ndReg_get => get
  use nuclearDatabase_inter,      only : nuclearDatabase

  implicit none
  private

  !!
  !! Transport operator that moves a particle with delta tracking
  !!
  type, public, extends(transportOperator) :: transportOperatorDT
  contains
    procedure :: transit => deltaTracking
    ! Override procedure
    procedure :: init

  end type transportOperatorDT

contains

  !!
  !! Performs delta tracking until a real collision point is found
  !!
  subroutine deltaTracking(self, p, tally, thisCycle, nextCycle)
    class(transportOperatorDT), intent(inout) :: self
    class(particle), intent(inout)            :: p
    type(tallyAdmin), intent(inout)           :: tally
    class(particleDungeon), intent(inout)     :: thisCycle
    class(particleDungeon), intent(inout)     :: nextCycle
    real(defReal)                             :: majorant_inv, sigmaT, distance
    character(100), parameter :: Here = 'deltaTracking (transportOperatorDT_class.f90)'

    ! Get majorant XS inverse: 1/Sigma_majorant
    majorant_inv = ONE / self % xsData % getTrackingXS(p, p % matIdx(), MAJORANT_XS)

    ! Should never happen! Prevents Inf distances
    if (abs(majorant_inv) > huge(majorant_inv)) call fatalError(Here, "Majorant is 0")

    DTLoop:do
      distance = -log( p% pRNG % get() ) * majorant_inv

      ! Move partice in the geometry
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
          print*, 'Particle location: ', p % rGlobal()
          call fatalError(Here, "Particle is in undefined material")

        ! Give error if the particle is in a region with overlapping cells
        case(OVERLAP_MAT)
          print*, 'Particle location: ', p % rGlobal()
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
  !! Initialise DT transport operator
  !!
  !! See transportOperator_inter for more details
  !!
  subroutine init(self, dict)
    class(transportOperatorDT), intent(inout) :: self
    class(dictionary), intent(in)             :: dict

    ! Initialise superclass
    call init_super(self, dict)

  end subroutine init


end module transportOperatorDT_class
