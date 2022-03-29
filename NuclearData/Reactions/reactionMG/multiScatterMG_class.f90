module multiScatterMG_class

  use numPrecision
  use endfConstants
  use genericProcedures,     only : fatalError, numToChar
  use RNG_class,             only : RNG
  use reactionHandle_inter,  only : reactionHandle
  use reactionMG_inter,      only : reactionMG
  use dataDeck_inter,        only : dataDeck
  use dictDeck_class,        only : dictDeck
  use dictionary_class,      only : dictionary

  implicit none
  private

  !!
  !! Public pointer cast
  !!
  public :: multiScatterMG_TptrCast
  public :: multiScatterMG_CptrCast

  !!
  !! Procedures to extend in subclasses
  !!
  public :: kill
  public :: buildFromDict

  !!
  !! Isotropic multiplicative scattering
  !!
  !! Simple object that contains data related to multiplicative scattering
  !! Multiplicative scattering is a G -> G scattering described with Legendre Polynomials
  !! and possible weight gain or reduction during scattering events (given by release)
  !!
  !! This implementation assumes no delay in 2nd-ary emission
  !!
  !! Public members:
  !!   scatterXSs -> vector of total P0 scattering XSs
  !!   P0         -> P0 scattering matrix [G_out; G_in]
  !!   prod       -> Production matrix    [G_out; G_in]
  !!
  !! Interface:
  !!   reactionMG interface
  !!   buildFromDict -> build the object from SCONE dictionary
  !!   sampleGout    -> Sample outgoin energy group given incident group
  !!   scatterXS     -> returns total P0 scattering XSs from any group
  !!   production    -> return production value for scattering from group G_in -> G_out
  !!
  !! NOTE:
  !!   MG data is presented as public members but it is safeter to access it through
  !!   access procedures. Break of encapsulation is mostly for conveniance when extending this
  !!   class, but note that it is bad practice to relay too much on implementation details of the
  !!   superclass!
  !!
  type, public, extends(reactionMG) :: multiScatterMG
    real(defReal),dimension(:),allocatable    :: scatterXSs
    real(defReal),dimension(:,:),allocatable  :: P0
    real(defReal), dimension(:,:),allocatable :: prod
  contains
    ! Superclass procedures
    procedure :: init
    procedure :: kill
    procedure :: release
    procedure :: releasePrompt
    procedure :: releaseDelayed
    procedure :: sampleDelayRate
    procedure :: sampleOut

    ! Local procedures
    procedure :: buildFromDict
    procedure, non_overridable :: sampleGout
    procedure, non_overridable :: scatterXS
    procedure, non_overridable :: production

  end type multiScatterMG

contains

  !!
  !! Initialise multiScatterMG from data deck
  !!
  !! See reactionHandle for details.
  !!
  !! Errors:
  !!   FatalError if given unsupported dataDeck.
  !!
  !! NOTE:
  !!   Ignores MT for now.
  !!
  subroutine init(self, data, MT)
    class(multiScatterMG), intent(inout) :: self
    class(dataDeck), intent(inout)       :: data
    integer(shortInt), intent(in)        :: MT
    character(100), parameter :: Here = 'init (multiScatterMG_class.f90)'

    ! Select dynamic type of data deck and build
    select type(data)
      type is(dictDeck)
        call self % buildFromDict(data % dict)

      class default
        call fatalError(Here,'multiScatterMG cannot be build from: '//data % myType())

    end select
  end subroutine init

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(multiScatterMG), intent(inout) :: self

    ! Clean memory
    if(allocated(self % scatterXSs)) deallocate(self % scatterXSs)
    if(allocated(self % P0))         deallocate(self % P0)
    if(allocated(self % prod))       deallocate(self % prod)

  end subroutine kill

  !!
  !! Returns number of particles produced on average by the scatter
  !!
  !! See reactionMG documentation for details
  !!
  !! NOTE:
  !!   Is not optimised for speed. It is possible that will be hardly used. Implement
  !!   but consider removing at the later date
  !!
  pure function release(self, G) result(N)
    class(multiScatterMG), intent(in) :: self
    integer(shortInt), intent(in)     :: G
    real(defReal)                     :: N

    if (G > 0 .and. G <= size(self % scatterXSs)) then
      ! Calculate XS weighted average of production
      N = sum(self % prod(:,G) * self % P0(:,G)) / self % scatterXSs(G)

    else
      N = ZERO

    end if

  end function release

  !!
  !! Returns number of particles produced instantly on average by the fission
  !!
  !! See reactionMG documentation for details
  !!
  pure function releasePrompt(self, G) result(N)
    class(multiScatterMG), intent(in)  :: self
    integer(shortInt), intent(in)      :: G
    real(defReal)                      :: N

    N = self % release(G)

  end function releasePrompt

  !!
  !! Returns number of particles produced with delay on average by the fission
  !!
  !! See reactionMG documentation for details
  !!
  !! NOTE:
  !!   For now returns always 0.0 !
  !!
  pure function releaseDelayed(self, G) result(N)
    class(multiScatterMG), intent(in)  :: self
    integer(shortInt), intent(in) :: G
    real(defReal)                 :: N

    N = ZERO

  end function releaseDelayed

  !!
  !! Sample the delay rate for the delayed particle
  !!
  !! See reactionMG documentation for details
  !!
  !! NOTE:
  !!   For now always returns 0.0 !
  !!
  function sampleDelayRate(self, G, rand) result(lambda)
    class(multiScatterMG), intent(in) :: self
    integer(shortInt), intent(in)     :: G
    class(RNG), intent(inout)         :: rand
    real(defReal)                     :: lambda

    lambda = ZERO

  end function sampleDelayRate

  !!
  !! Sample outgoing particle
  !!
  !! See reactionMG documentation for details
  !!
  subroutine sampleOut(self, mu, phi, G_out, G_in, rand)
    class(multiScatterMG), intent(in)   :: self
    real(defReal), intent(out)     :: mu
    real(defReal), intent(out)     :: phi
    integer(shortInt), intent(out) :: G_out
    integer(shortInt), intent(in)  :: G_in
    class(RNG), intent(inout)      :: rand
    character(100),parameter :: Here = 'sampleOut (multiScatterMG_class.f90)'

    ! Sample G_out
    G_out = self % sampleGout(G_in, rand)

    ! Sample deflection
    mu  = TWO * rand % get() - ONE
    phi = TWO_PI * rand % get()

  end subroutine sampleOut

  !!
  !! Sample outgoing energy group
  !!
  !! Args:
  !!   G_in [in]   -> incident energy group
  !!   rand [inout]-> random number generator
  !!
  !! Result:
  !!   Outgoing energy group G_out
  !!
  !! Errors:
  !!   fatalError if G_in is out of bounds
  !!
  function sampleGout(self, G_in, rand) result(G_out)
    class(multiScatterMG), intent(in) :: self
    integer(shortInt), intent(in)     :: G_in
    class(RNG), intent(inout)         :: rand
    integer(shortInt)                 :: G_out
    real(defReal)                     :: rem
    character(100),parameter :: Here = 'sampleGout (multiScatterMG_class.f90)'

    ! Check range
    if (G_in < 0 .or. G_in > size(self % scatterXSs)) then
      call fatalError(Here, 'Invalid incident group number: '//numToChar(G_in))
    end if

    ! Perform sampling
    rem = rand % get() * self % scatterXSs(G_in)
    do G_out = 1,size(self % scatterXSs)
      rem = rem - self % P0(G_out, G_in)
      if(rem < ZERO) return
    end do

    ! Sampling failed
    call fatalError(Here,'WTF? Sampling failed. Wrong scatter XS or random number above 1?')

  end function sampleGout

  !!
  !! Return total P0 scattering Cross-Section for group G_in
  !!
  !! Args:
  !!   G [in] -> incident energy group
  !!
  !! Result:
  !!   Sum of P0 scattering XSs over all outgoing groups
  !!
  !! Error:
  !!   For G_in out of bounds returns 0.0
  !!
  elemental function scatterXS(self, G) result(xs)
    class(multiScatterMG), intent(in) :: self
    integer(shortInt), intent(in)     :: G
    real(defReal)                     :: xs

    if( G < 0 .or. G > size(self % scatterXSs)) then
      xs = ZERO
    else
      xs = self % scatterXSs(G)
    end if

  end function scatterXS

  !!
  !! Return average production for scattering from G_in to G_out
  !!
  !! Args:
  !!   G_in [in] -> Incident energy group
  !!   G_out  [in]  -> Outgoing energy group
  !!
  !! Result:
  !!   Average particle production for transition from G_in to G_out
  !!
  !! Error:
  !!   If either G_out or G_in is out of bounds returns 0.0
  !!
  elemental function production(self, G_in, G_out) result (prod)
    class(multiScatterMG), intent(in) :: self
    integer(shortInt), intent(in)     :: G_in
    integer(shortInt), intent(in)     :: G_out
    real(defReal)                     :: prod

    ! Check range
    if( min(G_out, G_in) > 0 .and. max(G_out,G_in) <=  size(self % scatterXSs)) then
      prod = self % prod(G_out, G_in)

    else
      prod = ZERO
    end if

  end function production

  !!
  !! Builds multiScatterMG from SCONE dictionary
  !!
  !! Args:
  !!   dict [in] -> dictionary that contains data
  !!
  !! Errors:
  !!   FatalError if number of groups is not +ve
  !!   FatalError if P1 or scatteringMultiplicity does not match # of groups
  !!
  subroutine buildFromDict(self, dict)
    class(multiScatterMG), intent(inout)    :: self
    class(dictionary), intent(in)           :: dict
    real(defReal),dimension(:),allocatable  :: temp
    integer(shortInt)                       :: nG
    character(100),parameter :: Here = 'buildFromDict (multiScatterMG_class.f90)'

    ! Read number of groups
    call dict % get(nG,'numberOfGroups')
    if(nG <= 0) call fatalError(Here, 'Not +ve number of energy groups')

    ! Read scattering matrix
    call dict % get(temp, 'P0')
    if( size(temp) /= nG*nG) then
      call fatalError(Here,'Invalid size of P0. Expected: '//numToChar(nG**2)//&
                           ' got: '//numToChar(size(temp)))
    end if
    self % P0 = reshape(temp,[nG, nG])

    ! Read production matrix
    call dict % get(temp, 'scatteringMultiplicity')
    if( size(temp) /= nG*nG) then
      call fatalError(Here,'Invalid size of scatteringMultiplicity. Expected: '//numToChar(nG**2)//&
                           ' got: '//numToChar(size(temp)))
    end if

    self % prod = reshape(temp,[nG, nG])

    ! Calculate P0 total scattering XSs
    ! Behold the GLORY of Fortran you lowly C++ slaves!
    self % scatterXSs = sum(self % P0, 1)

  end subroutine buildFromDict

  !!
  !! Cast reactionHandle pointer to multiScatterMG pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class reactionHandle
  !!
  !! Result:
  !!   Null is source is not of multiScatterMG type
  !!   Target points to source if source is multiScatterMG type
  !!
  pure function multiScatterMG_TptrCast(source) result(ptr)
    class(reactionHandle), pointer, intent(in) :: source
    type(multiScatterMG), pointer              :: ptr

    select type(source)
      type is(multiScatterMG)
        ptr => source

      class default
        ptr => null()
    end select

  end function multiScatterMG_TptrCast

  !!
  !! Cast reactionHandle pointer to multiScatterMG Class pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class reactionHandle
  !!
  !! Result:
  !!   Null is source is not of multiScatterMG class
  !!   Target points to source if source is multiScatterMG type
  !!
  pure function multiScatterMG_CptrCast(source) result(ptr)
    class(reactionHandle), pointer, intent(in) :: source
    class(multiScatterMG), pointer              :: ptr

    select type(source)
      class is(multiScatterMG)
        ptr => source

      class default
        ptr => null()
    end select

  end function multiScatterMG_CptrCast

end module multiScatterMG_class
