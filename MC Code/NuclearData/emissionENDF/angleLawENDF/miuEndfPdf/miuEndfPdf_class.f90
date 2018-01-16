module miuEndfPdf_class

  use numPrecision
  use RNG_class, only : RNG

  implicit none
  private

  type,abstract, public :: miuEndfPdf
    !! Abstract class for all polymorphic objects containing probability density function (pdf)
    !! for secondary emission angle cosine miu.
    private
  contains
    procedure(sample), deferred        :: sample
    procedure(probabilityOf), deferred :: probabilityOf
  end type miuEndfPdf

  abstract interface

    function sample(self,rand) result(miu)
      import :: RNG,&
                miuEndfPdf, &
                defReal
      class(miuEndfPdf),intent(in) :: self
      class(RNG),intent(inout)     :: rand
      real(defReal)                :: miu
    end function sample

    function probabilityOf(self,miu) result(probability)
      import :: defReal,&
                miuEndfPdf
      class(miuEndfPdf), intent(in) :: self
      real(defReal), intent(in)     :: miu
      real(defReal)                 :: probability
    end function probabilityOf
  end interface

  type,public :: miuEndfPdf_ptr
    !! Pointer wrapper object on miuEndfPdf so arrays of pointers can be created
      private
      class(miuEndfPdf),pointer :: ptr => null()
    contains
      generic :: assignment(=) => assignPointer_ptr, assignPointer
      procedure :: sample => sample_ptr
      procedure :: probabilityOf => probabilityOf_ptr

      procedure, private :: assignPointer_ptr
      procedure, private :: assignPointer
  end type miuEndfPdf_ptr

contains
  subroutine assignPointer_ptr(LHS,RHS)
    class(miuEndfPdf_ptr),intent(out) :: LHS
    type(miuEndfPdf_ptr),intent(in)   :: RHS

    if(associated(LHS % ptr)) deallocate(LHS % ptr)
    LHS % ptr => RHS % ptr

  end subroutine

  subroutine assignPointer(LHS,RHS)
    class(miuEndfPdf_ptr),intent(out)    :: LHS
    class(miuEndfPdf),pointer,intent(in) :: RHS

    if(associated(LHS % ptr)) deallocate(LHS % ptr)
    LHS % ptr => RHS
  end subroutine assignPointer

  function sample_ptr(self,rand) result (miu)
    class(miuEndfPdf_ptr),intent(in) :: self
    class(RNG),intent(inout)         :: rand
    real(defReal)                    :: miu

    miu = self % ptr % sample(rand)
  end function sample_ptr

  function probabilityOf_ptr(self,miu) result(probability)
    class(miuEndfPdf_ptr),intent(in) :: self
    real(defReal),intent(in)         :: miu
    real(defReal)                    :: probability

    probability = self % ptr % probabilityOf(miu)

  end function probabilityOf_ptr


end module miuEndfPdf_class
