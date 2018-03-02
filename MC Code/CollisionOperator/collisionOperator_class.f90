module collisionOperator_class

  use numPrecision
  use endfConstants
  use genericProcedures, only : fatalError
  use RNG_class,         only : RNG
  use particle_class,    only : particle
  use byNucNoMT_class,   only : byNucNoMT

  ! Cross-section packages to interface with nuclear data
  use matNucCDF_class,   only : matNucCDF
  use xsMainCDF_class,   only : xsMainCDF

  implicit none
  private

  type, public :: collisionOperator
    private
    class(byNucNoMT), pointer :: xsData => null()
    class(RNG), pointer       :: locRNG => null()
  contains
    procedure :: attachXsData
    procedure :: collide
  end type collisionOperator

contains

  !!
  !! Initialise XS operator by providing pointer to XS data block
  !!
  subroutine attachXsData(self,xsData)
    class(collisionOperator), intent(inout) :: self
    class(byNucNoMT),pointer, intent(in)    :: xsData
    character(100), parameter               :: Here =' attachXsData (collisionOperator_class.f90)'

    if(.not.associated(xsData)) call fatalError(Here,'Allocated xs data must be provided')

    self % xsData => xsData

  end subroutine attachXSData


  !!
  !! Subroutine to collide a neutron. Chooses collision nuclide and main reaction channel.
  !! Calls approperiate procedure to change neutron state
  !!
  subroutine collide(self,p)
    class(collisionOperator), intent(inout) :: self
    type(particle), intent(inout)           :: p
    real(defReal)                           :: E
    integer(shortInt)                       :: matIdx
    integer(shortInt)                       :: nucIdx
    integer(shortInt)                       :: MT
    type(matNucCDF),pointer                 :: nuclideCDF
    type(xsMainCDF),pointer                 :: reactionCDF
    real(defReal)                           :: r

    ! Retrive Pointer to the random number generator
    self % locRNG => p % pRNG

    ! Load neutron energy and material
    E = p % E
    matIdx = p % matIdx

    ! Select collision nuclide
    call self % xsData % getMatNucCDF(nuclideCDF, E, matIdx)

    r = self % locRNG % get()
    nucIdx = nuclideCDF % invert(r)

    ! Select Main reaction channel
    call self % xsData % getMainNucCdf(reactionCDF, E, nucIdx)

    r = self % locRNG % get()
    MT = reactionCDF % invert(r)

    ! Call procedure to do reaction processing
    select case (MT)
      case(anyScatter)
     !   call self % performScattering(p,nucIdx)

      case(anyCapture)
     !   call self % performCapture(p,nucIdx)

      case(anyFission)
     !   call self % performFission(p,nucIdx)

    end select

  end subroutine collide

  !!
  !! Change particle state in scattering reaction
  !!
  subroutine performScattering(self,p,nucIdx)
    class(collisionOperator), intent(inout) :: self
    type(particle), intent(inout)           :: p
    integer(shortInt),intent(in)            :: nucIdx

    ! Select if target is stationary or not

     ! Call prodecure for stationary treatment

     ! Call procedure fo moving treatment



  end subroutine performScattering

  subroutine scatterFromStationary(self,p,E,nucIdx)
    class(collisionOperator), intent(inout) :: self
    type(particle), intent(inout)           :: p
    real(defReal), intent(in)               :: E       ! Neutron energy
    integer(shortInt),intent(in)            :: nucIdx  ! Target nuclide index

    ! Find coordinate frame
    ! Sample mu and outgoing energy
    ! Change mu and outgoing energy to LAB frame if given in CM
    ! Check energy cutoff and perhaps kill neutron
    ! Sample azimuthal deflection
    ! Turn particle

  end subroutine scatterFromStationary

  subroutine scatterFromMoving(self,p,E,nucIdx)
    class(collisionOperator), intent(inout) :: self
    type(particle), intent(inout)           :: p
    real(defReal), intent(in)               :: E       ! Neutron energy
    integer(shortInt),intent(in)            :: nucIdx  ! Target nuclide index

    ! Sample velocity of target
    ! Change coordinate to TARGET frame
    ! scatter from stationary
    ! Change back to LAB frame


  end subroutine scatterFromMoving




  !!
  !! Change particle state in capture reaction
  !!
  subroutine performCapture(self,p,nucIdx)
    class(collisionOperator), intent(inout) :: self
    type(particle), intent(inout)           :: p
    integer(shortInt),intent(in)            :: nucIdx

    p % isDead =.true.

  end subroutine performCapture

  !!
  !! Change particle state in fission reaction
  !!
  subroutine performFission(self,p,nucIdx)
    class(collisionOperator), intent(inout) :: self
    type(particle), intent(inout)           :: p
    integer(shortInt),intent(in)            :: nucIdx

    p % isDead =.true.

  end subroutine performFission



end module collisionOperator_class
