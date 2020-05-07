module muEndfPdfSlot_class

  use numPrecision
  use genericProcedures, only : fatalError
  use RNG_class,         only : RNG
  use aceCard_class,     only : aceCard
  use muEndfPdf_inter,   only : muEndfPdf

  ! Implementations
  use isotropicMu_class, only : isotropicMu
  use equiBin32Mu_class, only : equiBin32Mu
  use tabularMu_class,   only : tabularMu

  implicit none
  private

  !!
  !! Slot to store polymorphic instances of muEndfPdfSlot
  !! This breaks rules of standard inter,slot,Factory SCONE pattern
  !! because factory is included insied the slot (init subroutine)
  !!
  type, public,extends(muEndfPdf) :: muEndfPdfSlot
    private
    class(muEndfPdf),allocatable :: slot
  contains
    ! Duplicate interface of the superclass
    procedure :: sample
    procedure :: probabilityOf
    procedure :: kill

    ! Define assignment
    procedure :: init
    procedure :: moveAllocFrom

  end type muEndfPdfSlot

contains

  !!
  !! Sample mu given random number generator
  !!
  function sample(self,rand) result(mu)
    class(muEndfPdfSlot),intent(in)  :: self
    class(RNG),intent(inout)         :: rand
    real(defReal)                    :: mu

    mu = self % slot % sample(rand)

  end function sample

  !!
  !! Give probability density of mu
  !! Does not check if mu is in <-1;1>
  !!
  function probabilityOf(self,mu) result(prob)
    class(muEndfPdfSlot), intent(in)  :: self
    real(defReal), intent(in)         :: mu
    real(defReal)                     :: prob

    prob = self % slot % probabilityOf(mu)

  end function probabilityOf

  !!
  !! Initialise slot to a specific type determined by type
  !!
  subroutine init(self, ACE, type)
    class(muEndfPdfSlot), intent(inout) :: self
    class(aceCard), intent(inout)       :: ACE
    character(*), intent(in)            :: type
    character(100), parameter :: Here ='init (muEndfPdfSlot_class.f90)'

    select case(type)
      case('isotropicMu')
        allocate(self % slot, source = isotropicMu())

      case('equiBin32Mu')
        allocate(self % slot, source = equiBin32Mu(ACE))

      case('tabularMu')
        allocate( self % slot, source = tabularMu(ACE))

      case default
        print '(A)','Available implementations:','isotropicMu', 'equiBin32Mu', 'tabularMu'
        call fatalError(Here, 'Unrecognised type of muEndfPdf requested: ' // trim(type) )

    end select


  end subroutine init

  !!
  !! Move allocation from RHS to LHS slot
  !!
  subroutine moveAllocFrom(LHS,RHS)
    class(muEndfPdfSlot), intent(inout) :: LHS
    type(muEndfPdfSlot), intent(inout)  :: RHS

    if(allocated(LHS % slot)) deallocate (LHS % slot)

    call move_alloc(RHS % slot, LHS % slot)

  end subroutine moveAllocFrom

  !!
  !! Return to uninitialised state
  !!
  elemental subroutine kill(self)
    class(muEndfPdfSlot), intent(inout) :: self

    ! Kill slot if allocated
    if(allocated(self % slot)) then
      call self % slot % kill()
      deallocate(self % slot)
    end if

  end subroutine kill

end module muEndfPdfSlot_class
