module multiScatterP1MG_class

  use numPrecision
  use endfConstants
  use genericProcedures,     only : fatalError, numToChar
  use legendrePoly_func,     only : sampleLegendre
  use RNG_class,             only : RNG
  use reactionHandle_inter,  only : reactionHandle
  use dataDeck_inter,        only : dataDeck
  use dictDeck_class,        only : dictDeck
  use dictionary_class,      only : dictionary
  use multiScatterMG_class,  only : multiScatterMG, kill_super => kill, &
                                    buildFromDict_super => buildFromDict



  implicit none
  private

  !!
  !! Public pointer cast
  !!
  public :: multiScatterP1MG_TptrCast

  !!
  !! MG multiplicative scattering with P1 information
  !!
  !! Public Members:
  !!   P1 -> Matrix with normalised P1 scattering coefficients (P1 /P0)
  !!
  !! Interface
  !!   multiScatterMG interface
  !!
  !! NOTE:
  !!   It is sufficient to extend concreate build procedure. We can still use init in the
  !!   superclass!
  !!
  type, public, extends(multiScatterMG) :: multiScatterP1MG
    real(defReal), dimension(:,:), allocatable :: P1

  contains
    ! Override some procedures
    procedure :: kill
    procedure :: sampleOut
    procedure :: buildFromDict
  end type multiScatterP1MG

contains

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(multiScatterP1MG), intent(inout) :: self

    ! Call superclass procedure
    call kill_super(self)

    ! Clean own memory
    if(allocated(self % P1)) deallocate(self % P1)

  end subroutine kill

  !!
  !! Sample outgoing particle
  !!
  !! See reactionMG documentation for details
  !!
  subroutine sampleOut(self, mu, phi, G_out, G_in, rand)
    class(multiScatterP1MG), intent(in)   :: self
    real(defReal), intent(out)     :: mu
    real(defReal), intent(out)     :: phi
    integer(shortInt), intent(out) :: G_out
    integer(shortInt), intent(in)  :: G_in
    class(RNG), intent(inout)      :: rand
    character(100),parameter :: Here = 'sampleOut (multiScatterMG_class.f90)'

    ! Sample G_out
    G_out = self % sampleGout(G_in, rand)

    ! Sample deflection
    mu  = sampleLegendre(self % P1(G_out, G_in), rand)
    phi = TWO_PI * rand % get()

  end subroutine sampleOut

  !!
  !! Builds multiScatterP1MG from SCONE dictionary
  !!
  !! Extends multiScatterMG procedure!
  !! See its documentation for extra details.
  !!
  !! Errors:
  !!   FatalError if size of P1 scattering matrix does not match numer of group
  !!
  subroutine buildFromDict(self, dict)
    class(multiScatterP1MG), intent(inout)  :: self
    class(dictionary), intent(in)           :: dict
    real(defReal),dimension(:),allocatable  :: temp
    integer(shortInt)                       :: nG
    character(100),parameter :: Here = 'buildFromDict (multiScatterMG_class.f90)'

    ! Call superclass procedure
    call buildFromDict_super(self, dict)

    ! Re-read number of groups
    call dict % get(nG,'numberOfGroups')

    ! Read P1 scattering matrix
    call dict % get(temp, 'P1')
    if( size(temp) /= nG*nG) then
      call fatalError(Here,'Invalid size of P1. Expected: '//numToChar(nG**2)//&
                           ' got: '//numToChar(size(temp)))
    end if
    self % P1 = reshape(temp,[nG, nG])

    ! Normalise P1 coefficients
    where (self % P0 /= ZERO)
      self % P1 = self % P1 / self % P0

    elsewhere
      self % P1 = ZERO

    end where

  end subroutine buildFromDict

  !!
  !! Cast reactionHandle pointer to multiScatterP1MG pointer
  !!
  !! Args:
  !!   source [in]    -> source pointer of class reactionHandle
  !!
  !! Result:
  !!   Null is source is not of multiScatterP1MG type
  !!   Target points to source if source is multiScatterP1MG type
  !!
  pure function multiScatterP1MG_TptrCast(source) result(ptr)
    class(reactionHandle), pointer, intent(in) :: source
    type(multiScatterP1MG), pointer              :: ptr

    select type(source)
      type is(multiScatterP1MG)
        ptr => source

      class default
        ptr => null()
    end select

  end function multiScatterP1MG_TptrCast


end module multiScatterP1MG_class
