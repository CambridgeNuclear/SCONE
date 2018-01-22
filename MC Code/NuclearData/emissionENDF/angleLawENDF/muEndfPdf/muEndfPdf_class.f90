module muEndfPdf_class

  use numPrecision
  use RNG_class, only : RNG

  implicit none
  private

  type,abstract, public :: muEndfPdf
    !! Abstract class for all polymorphic objects containing probability density function (pdf)
    !! for secondary emission angle cosine mu.
    private
  contains
    procedure(sample), deferred        :: sample
    procedure(probabilityOf), deferred :: probabilityOf
  end type muEndfPdf

  abstract interface

    function sample(self,rand) result(mu)
      import :: RNG,&
                muEndfPdf, &
                defReal
      class(muEndfPdf),intent(in)  :: self
      class(RNG),intent(inout)     :: rand
      real(defReal)                :: mu
    end function sample

    function probabilityOf(self,mu) result(probability)
      import :: defReal,&
                muEndfPdf
      class(muEndfPdf), intent(in)  :: self
      real(defReal), intent(in)     :: mu
      real(defReal)                 :: probability
    end function probabilityOf
  end interface

  type,public :: muEndfPdf_ptr
    !! Pointer wrapper object on muEndfPdf so arrays of pointers can be created
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
  subroutine assignPointer_ptr(LHS,RHS)
    class(muEndfPdf_ptr),intent(out) :: LHS
    type(muEndfPdf_ptr),intent(in)   :: RHS

    if(associated(LHS % ptr)) deallocate(LHS % ptr)
    LHS % ptr => RHS % ptr

  end subroutine

  subroutine assignPointer(LHS,RHS)
    class(muEndfPdf_ptr),intent(out)    :: LHS
    class(muEndfPdf),pointer,intent(in) :: RHS

    if(associated(LHS % ptr)) deallocate(LHS % ptr)
    LHS % ptr => RHS
  end subroutine assignPointer

  function sample_ptr(self,rand) result (mu)
    class(muEndfPdf_ptr),intent(in)  :: self
    class(RNG),intent(inout)         :: rand
    real(defReal)                    :: mu

    mu = self % ptr % sample(rand)
  end function sample_ptr

  function probabilityOf_ptr(self,mu) result(probability)
    class(muEndfPdf_ptr),intent(in)  :: self
    real(defReal),intent(in)         :: mu
    real(defReal)                    :: probability

    probability = self % ptr % probabilityOf(mu)

  end function probabilityOf_ptr


end module muEndfPdf_class
