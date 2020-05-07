!!
!! This module contains a number of classes that can contain results of a tallyClerk
!!   These results are for the purpose of interaction with Physics Packages (not output)
!!
!!   tallyResult       is an abstract class to create family of result types
!!   tallyResultScalar is a simple result consisting of a single scalar
!!   tallyResultEmpty  is a null object returned by clerks where definition of tallyResult was
!!                     not implemented
!!
!! More complex results, should be written alongside tallyClerks that produce them.
!! There is no general N-dimensional result type, becouse it would be an object of very high
!! complexity and it is safe to assume that a Physics Package interacting with the result
!! "knows" what kind of result to expect, so a specialised tallyResult for a particular output
!! will perform much better in terms of clarity
!!
module tallyResult_class

  use numPrecision

  implicit none
  private

  !!
  !! Empty abstract class to allow polymorphism of all result subclasses
  !!
  type, public,abstract :: tallyResult

  end type tallyResult

  !!
  !! Very simple class for a single scalar result. Can be used by multiple
  !! Use default constructor e.g.:
  !!
  !! var = tallyResultScalar(name, value, STD)
  !!
  type,public, extends(tallyResult) :: tallyResultScalar
    character(nameLen) :: clerkName
    real(defReal)      :: value
    real(defReal)      :: STD
  end type tallyResultScalar

  !!
  !! Class that is returned when there is no defined result
  !!
  type,public,extends(tallyResult) :: tallyResultEmpty

  end type tallyResultEmpty

end module tallyResult_class
