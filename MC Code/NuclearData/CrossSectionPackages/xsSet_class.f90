module xsSet_class

  use numPrecision

  implicit none
  private
  !!
  !! Type to store a reaction cross-sections. It should be possible to acess xss by their
  !! associated MT number. However, child classes will break encapsulation not to overcomplicate
  !! the interfaces. There is no point to hide a lumped scattering or capture cross section behind
  !! MT number. As the result, this abstract class mey quite unnecessary. It may be removed at the
  !! later stage.
  !!
  type,abstract, public :: xsSet
    private
  contains
    procedure(xsOf), deferred :: xsOf
  end type xsSet


  abstract interface

    function xsOf(self,MT) result(xs)
      import :: xsSet, &
                shortInt, &
                defReal
      class(xsSet), intent(in)      :: self
      integer(shortInt), intent(in) :: MT
      real(defReal)                 :: xs

    end function xsOf

  end interface


contains

    
end module xsSet_class
