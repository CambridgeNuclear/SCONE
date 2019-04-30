module multiMap_test
  use numPrecision
  use pFUnit_mod
  use particle_class,          only : particleState
  use dictionary_class,        only : dictionary
  use multiMap_class,          only : multiMap

  implicit none


@testCase
  type, extends(TestCase) :: test_multiMap
    type(multiMap) :: map

  contains
    procedure :: setUp
    procedure :: tearDown
  end type test_multiMap

contains

  !!
  !! Sets up test_intMap object we can use in a number of tests
  !!
  subroutine setUp(this)
    class(test_multiMap), intent(inout) :: this
    type(dictionary)                    :: tempDict


  end subroutine setUp

  !!
  !! Kills test_intMap object we can use in a number of tests
  !!
  subroutine tearDown(this)
    class(test_multiMap), intent(inout) :: this


  end subroutine tearDown

!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! PROPER TESTS BEGIN HERE
!!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  !!
  !! Test that maps performs as expected
  !!
@Test
  subroutine multiMapping(this)
    class(test_multiMap), intent(inout) :: this
    type(particleState) :: state


  end subroutine multiMapping

end module multiMap_test
