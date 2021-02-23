module universe_test
  use numPrecision
  use charMap_class,  only : charMap
  use universe_inter, only : charToFill
  use pFUnit_mod
  implicit none


  !!
  !! Note that universe is abstract thus it canot be tested by itself
  !!
  !! Tests for universe non-overridable procedures are in cellUniverse_test
  !!

contains

  !!
  !! Test charToFill
  !!
  !! Since universe is abstract it cannot be tested by itself (only via its
  !! subclasses)
  !!
@Test
  subroutine test_charToFill()
    type(charMap)      :: mats
    character(nameLen) :: name
    character(100), parameter :: Here = 'parentScope'

    ! Load some material names and their (fake) matIdxs
    name = 'mat13'
    call mats % add(name, 13)

    name = 'city17'
    call mats % add(name, 17)

    name = 'mat47'
    call mats % add(name, 47)

    ! Test material converstion
    name = 'mat13'
    @assertEqual(13, charToFill(name, mats, Here))

    name = 'city17'
    @assertEqual(17, charToFill(name, mats, Here))

    ! Test universe ID conversion
    ! By convention returns -ve uniID
    ! NO SPACES IN THE CHAR!!!
    name = 'u<87>'
    @assertEqual(-87, charToFill(name, mats, Here))

    name = 'u<133>'
    @assertEqual(-133, charToFill(name, mats, Here))

    name = 'u<5>'
    @assertEqual(-5, charToFill(name, mats, Here))

  end subroutine test_charToFill



end module universe_test
