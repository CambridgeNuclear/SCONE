module levelScattering_class

  use numPrecision
  use genericProcedures,   only : fatalError, numToChar
  use aceCard_class,       only : aceCard
  use RNG_class,           only : RNG
  use energyLawENDF_inter, only : energyLawENDF

  implicit none
  private

  interface levelScattering
    module procedure new_levelScattering
    module procedure new_levelScattering_fromACE
  end interface


  !!
  !! Level Scattering (ENDF Energy Law 3)
  !! Not a probability distribution.
  !! Outgoing energy is a function of incident energy
  !! Description in MCNP Manual Appendix F TABLE F-14 b)
  !!
  !! Private members:
  !!   LDAT1 -> (A+1)/A * (-Q) which is +ve for endothermic scattering and -ve for
  !!       exothermic scattering
  !!   LDAT2 -> (A/(A+1))^2 must be in [0,1]
  !!
  type, public,extends(energyLawENDF) :: levelScattering
    private
    real(defReal) :: LDAT1 = ZERO
    real(defReal) :: LDAT2 = ZERO ! (A/(A+1))^2
  contains
    ! Interface procedures
    procedure :: sample
    procedure :: probabilityOf
    procedure :: kill

    ! Instance procedures
    procedure :: init
  end type levelScattering

contains

  !!
  !! "Sample" outgoing energy given random number generator and incedent energy
  !! Not much to sample. This type of scattering is deterministic
  !!
  function sample(self,E_in,rand) result (E_out)
    class(levelScattering), intent(in) :: self
    real(defReal), intent(in)          :: E_in
    class(RNG), intent(inout)          :: rand
    real(defReal)                      :: E_out

    ! Calculate outgoing energy
    E_out = self % LDAT2 * (E_in- self % LDAT1)

  end function sample

  !!
  !! Given probability of outgoing energy given incedent energy
  !! Given that scattering is deterministic probability is a delta function
  !!
  function probabilityOf(self,E_out,E_in) result (prob)
    class(levelScattering), intent(in) :: self
    real(defReal), intent(in)          :: E_out,E_in
    real(defReal)                      :: prob
    real(defReal)                      :: E_temp

    ! Check if E_out is equal to outgoing energy for a given E_in
    E_temp = self % LDAT2 * (E_in- self % LDAT1)
    if(E_out == E_temp) then
      prob = ONE

    else
      prob = ZERO

    end if
  end function probabilityOf

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(levelScattering), intent(inout) :: self

    self % LDAT1 = ZERO
    self % LDAT2 = ZERO

  end subroutine kill

  !!
  !! Initialise level scattering
  !!
  !! Args:
  !!   LDAT1 [in] -> real number
  !!   LDAT2 [in] -> Value of (A/A+1)**2 in [0,1]
  !!
  !! Errors:
  !!   Returns fatalError if LDAT2 is outside [0,1]
  !!
  subroutine init(self,LDAT1,LDAT2)
    class(levelScattering), intent(inout) :: self
    real(defReal), intent(in)             :: LDAT1
    real(defReal), intent(in)             :: LDAT2
    character(100),parameter :: Here='init (levelScattering_class.f90)'

    ! Perform sanity checks
    if( LDAT2 <  ZERO ) call fatalError(Here,'LDAT2 is -ve:' // numToChar(LDAT2))
    if( LDAT2 >= ONE  ) call fatalError(HEre,'LDAT2 is >= 1.0:' // numToChar(LDAT2))

    ! Assign values
    self % LDAT1 = LDAT1
    self % LDAT2 = LDAT2

  end subroutine init

  !!
  !! Constructor
  !! See init for details
  !!
  function new_levelScattering(LDAT1,LDAT2) result(new)
    real(defReal), intent(in)             :: LDAT1
    real(defReal), intent(in)             :: LDAT2
    type(levelScattering)                 :: new

    ! Initialise
    call new % init(LDAT1, LDAT2)

  end function new_levelScattering

  !!
  !! Constructor from ACE
  !!
  !! Args:
  !!   ACE [inout] -> ACE card set to the beginning of data
  !!
  !! Errors:
  !!   Will return invalid value without error if ACE card is set to
  !!   wrong location.
  !!
  function new_levelScattering_fromACE(ACE) result(new)
    type(aceCard), intent(inout)          :: ACE
    type(levelScattering)                 :: new
    real(defReal)                         :: LDAT1, LDAT2

    ! Read data
    LDAT1 = ACE % readReal()
    LDAT2 = ACE % readReal()

    ! Initialise
    call new % init(LDAT1, LDAT2)

  end function new_levelScattering_fromACE


end module levelScattering_class
