module muEndfPdf_inter

  use numPrecision
  use RNG_class,      only : RNG

  implicit none
  private

  !!
  !! Abstract interface  for all polymorphic objects containing probability density function (pdf)
  !! for secondary emission angle cosine mu.
  !!
  type,abstract, public :: muEndfPdf
    private
  contains
    procedure(sample), deferred        :: sample
    procedure(probabilityOf), deferred :: probabilityOf
  end type muEndfPdf

  abstract interface

    !!
    !! Sample mu given random number generator
    !!
    function sample(self,rand) result(mu)
      import :: RNG,&
                muEndfPdf, &
                defReal
      class(muEndfPdf),intent(in)  :: self
      class(RNG),intent(inout)     :: rand
      real(defReal)                :: mu
    end function sample

    !!
    !! Give probability density of mu
    !! Does not check if mu is in <-1;1>
    !!
    function probabilityOf(self,mu) result(probability)
      import :: defReal,&
                muEndfPdf
      class(muEndfPdf), intent(in)  :: self
      real(defReal), intent(in)     :: mu
      real(defReal)                 :: probability
    end function probabilityOf
  end interface

  !! *** OBSOLETE WILL BE REPLACED WITH SLOTS SOON
  !! Pointer wrapper object on muEndfPdf so arrays of pointers can be created
  !!
  type,public :: muEndfPdf_ptr
      private
      class(muEndfPdf),pointer :: ptr => null()
    contains
      generic :: assignment(=) => assignPointer_ptr, assignPointer
      procedure :: sample => sample_ptr
      procedure :: probabilityOf => probabilityOf_ptr

      procedure, private :: assignPointer_ptr
      procedure, private :: assignPointer
  end type muEndfPdf_ptr

contains

  !!
  !! Overload assignment.
  !! Copy pointer from pointer wrapper to pointer wrapper
  !!
  subroutine assignPointer_ptr(LHS,RHS)
    class(muEndfPdf_ptr),intent(out) :: LHS
    type(muEndfPdf_ptr),intent(in)   :: RHS

    if(associated(LHS % ptr)) deallocate(LHS % ptr)
    LHS % ptr => RHS % ptr

  end subroutine

  !!
  !! Overload assignment.
  !! Copy pointer to pointer wrapper
  !!
  subroutine assignPointer(LHS,RHS)
    class(muEndfPdf_ptr),intent(out)    :: LHS
    class(muEndfPdf),pointer,intent(in) :: RHS

    if(associated(LHS % ptr)) deallocate(LHS % ptr)
    LHS % ptr => RHS
  end subroutine assignPointer

  !!
  !! Tranfer sample call to content
  !!
  function sample_ptr(self,rand) result (mu)
    class(muEndfPdf_ptr),intent(in)  :: self
    class(RNG),intent(inout)         :: rand
    real(defReal)                    :: mu

    mu = self % ptr % sample(rand)
  end function sample_ptr

  !!
  !! Tranfer probabilityOf call to content
  !!
  function probabilityOf_ptr(self,mu) result(probability)
    class(muEndfPdf_ptr),intent(in)  :: self
    real(defReal),intent(in)         :: mu
    real(defReal)                    :: probability

    probability = self % ptr % probabilityOf(mu)

  end function probabilityOf_ptr


end module muEndfPdf_inter
