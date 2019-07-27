module kalbachPdf_class

  use numPrecision
  use genericProcedures,  only : fatalError
  use aceCard_class,      only : aceCard
  use RNG_class,          only : RNG
  use kalbachTable_class, only : kalbachTable

  implicit none
  private

  interface kalbachPdf
    module procedure new_kalbachPdf_withPDF
    module procedure new_kalbachPdf_withCDF
    module procedure new_kalbachPdf_fromACE
  end interface

  !!
  !! Probability distribution table for kalbach data at a single energy point
  !!
  type, public :: kalbachPdf
    private
    type(kalbachTable) :: table
  contains
    ! Functional procedures
    procedure :: sample
    procedure :: bounds
    procedure :: probabilityOf
    procedure :: kill

    ! Initialisation procedures
    generic           :: init          => init_withPDF, init_withCDF
    procedure,private :: init_withPDF
    procedure,private :: init_withCDF
  end type kalbachPdf

contains

  !!
  !! Sample outgoing angle mu and energy E_out given  random number generator
  !!
  subroutine sample(self,mu,E_out,rand)
    class(kalbachPdf), intent(in)   :: self
    real(defReal),intent(out)       :: mu
    real(defReal),intent(out)       :: E_out
    class(RNG),intent(inout)        :: rand
    real(defReal)                   :: R,A,T
    real(defReal)                   :: r1,r2,r3

    ! Generate random number
    r1 = rand % get()

    ! Sample outgoing energy
    call self % table % sample(r1,E_out,R,A)

    ! Sample mu -> scheme copied from MCNP manual Chapter 2
    r2 = rand % get()
    r3 = rand % get()

    if( r2 >= R) then
      T = (TWO * r3 - ONE)*sinh(A)
      mu = log(T+sqrt(T*T+ONE))/A

    else
      mu = log(r3 * exp(A) + (ONE-r3) * exp(-A))/A

    end if

  end subroutine sample

  !!
  !! Return energy bounds of the probability distribution
  !!
  subroutine bounds(self,E_min,E_max)
    class(kalbachPdf), intent(in)    :: self
    real(defReal), intent(out)       :: E_min
    real(defReal), intent(out)       :: E_max

    call self % table % bounds(E_min, E_max)

  end subroutine bounds

  !!
  !! Returns probability that neutron was emitted at angle mu and energy E_out
  !!
  function probabilityOf(self,mu,E_out) result (prob)
    class(kalbachPdf), intent(in)   :: self
    real(defReal),intent(in)        :: mu
    real(defReal),intent(in)        :: E_out
    real(defReal)                   :: prob
    real(defReal)                   :: P_Eout, R, A, P_mu

    ! Obtain Kalbach parameters and probability of E_out
    call self % table % probabilityOf(E_out, P_Eout, R, A)

    ! Calculate probability of mu
    P_mu = HALF * A/sinh(A)*(cosh(A*mu) + R*sinh(A*mu))

    ! Combined probability
    prob =  P_mu * P_Eout

  end function probabilityOf

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(kalbachPdf), intent(inout) :: self

    call self % table % kill()

  end subroutine kill


  !!
  !! Initialise with PDF only
  !!
  subroutine init_withPDF(self,E,pdf,R,A,interFlag)
    class(kalbachPdf), intent(inout)       :: self
    real(defReal),dimension(:),intent(in)  :: E
    real(defReal),dimension(:),intent(in)  :: pdf
    real(defReal),dimension(:),intent(in)  :: R
    real(defReal),dimension(:),intent(in)  :: A
    integer(shortInt),intent(in)           :: interFlag
    character(100),parameter :: Here ='init_withPDF (kalbachPdf_class.f90)'

    ! Perform checks
    if(any( E < 0.0 ) ) call fatalError(Here,'E contains -ve values')
    
    ! Initialise table
    call self % table % init(E,pdf,R,A,interFlag)

  end subroutine init_withPDF

  !!
  !! Initialise with PDF and CDF
  !! Does NOT check if PDF and CDF are consistant
  !!
  subroutine init_withCDF(self,E,pdf,cdf,R,A,interFlag)
    class(kalbachPdf), intent(inout)       :: self
    real(defReal),dimension(:),intent(in)  :: E
    real(defReal),dimension(:),intent(in)  :: pdf
    real(defReal),dimension(:),intent(in)  :: cdf
    real(defReal),dimension(:),intent(in)  :: R
    real(defReal),dimension(:),intent(in)  :: A
    integer(shortInt),intent(in)           :: interFlag
    character(100),parameter :: Here ='init_withCDF (kalbachPdf_class.f90)'

    ! Perform checks
    if(any( E < 0.0 ) ) call fatalError(Here,'E contains -ve values')

    ! Initialise table
    call self % table % init(E,pdf,cdf,R,A,interFlag)

  end subroutine init_withCDF

  !!
  !! Constructor with PDF
  !!
  function new_kalbachPdf_withPDF(E,pdf,R,A,interFlag) result(new)
    real(defReal),dimension(:),intent(in)  :: E
    real(defReal),dimension(:),intent(in)  :: pdf
    real(defReal),dimension(:),intent(in)  :: R
    real(defReal),dimension(:),intent(in)  :: A
    integer(shortInt),intent(in)           :: interFlag
    type(kalbachPdf)                       :: new

    ! Initialise
    call new % init(E,pdf,R,A,interFlag)

  end function new_kalbachPdf_withPDF

  !!
  !! Constructor with PDF and CDF
  !!
  function new_kalbachPdf_withCDF(E,pdf,cdf,R,A,interFlag) result(new)
    real(defReal),dimension(:),intent(in)  :: E
    real(defReal),dimension(:),intent(in)  :: pdf
    real(defReal),dimension(:),intent(in)  :: cdf
    real(defReal),dimension(:),intent(in)  :: R
    real(defReal),dimension(:),intent(in)  :: A
    integer(shortInt),intent(in)           :: interFlag
    type(kalbachPdf)                       :: new

    ! Initialise
    call new % init(E,pdf,cdf,R,A,interFlag)

  end function new_kalbachPdf_withCDF

  !!
  !! Constructor from ACE
  !! aceCard read head needs to point to the beginning of data
  !!
  function new_kalbachPdf_fromACE(ACE) result(new)
    type(aceCard), intent(inout)            :: ACE
    type(kalbachPdf)                        :: new
    real(defReal),dimension(:),allocatable  :: E
    real(defReal),dimension(:),allocatable  :: pdf
    real(defReal),dimension(:),allocatable  :: cdf
    real(defReal),dimension(:),allocatable  :: R
    real(defReal),dimension(:),allocatable  :: A
    integer(shortInt)                       :: interFlag
    integer(shortInt)                       :: N

    ! Read interpolation flag
    interFlag = ACE % readInt()

    ! Read number of points in distribution
    N = ACE % readInt()

    ! Read rest of the data
    E   = ACE % readRealArray(N) ! Energy grid
    pdf = ACE % readRealArray(N) ! Probability Density Function
    cdf = ACE % readRealArray(N) ! Cumulative Distribution Function
    R   = ACE % readRealArray(N) ! Precompound fraction
    A   = ACE % readRealArray(N) ! Angular distribution slope

    ! initialise
    call new % init(E,pdf,cdf,R,A,interFlag)

  end function new_kalbachPdf_fromACE

end module kalbachPdf_class
