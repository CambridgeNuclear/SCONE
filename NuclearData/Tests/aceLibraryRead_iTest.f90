module aceLibraryRead_iTest

  use numPrecision
  use aceLibrary_mod, only : load, new_neutronACE
  use aceCard_class,  only : aceCard
  use pFUnit_mod

  implicit none

contains

  !!
  !! Test reading ACE Library
  !!   Make sure that a libery file is read without any problems and that all ACE cards
  !!   can be read as well
  !!
@Test
  subroutine testReadingACELibrary()
    type(aceCard)      :: ACE
    character(nameLen) :: ZAID
    real(defReal), parameter :: TOL = 1.0E-6


    ! Load library
    call load('./IntegrationTestFiles/testLib')

    ! Load U-233
    ZAID = '92233.03'
    call new_neutronACE(ACE,ZAID)

    ! Verify
    @assertEqual('92233.03c', ACE % ZAID)
    @assertEqual(231.038000_defReal, ACE % AW, TOL * 231.038000_defReal)
    @assertEqual(25, ACE % numMT())

    ! Load H-1 at 300 K
    ZAID = '1001.03'
    call new_neutronACE(ACE,ZAID)
    @assertEqual('1001.03c', ACE % ZAID)
    @assertEqual(0.999170_defReal, ACE % AW, TOL)
    @assertEqual(2.5852E-08 , ACE % TZ ,TOL * 2.5852E-08)
    @assertEqual(1, ACE % numMT())


  end subroutine testReadingACELibrary


end module aceLibraryRead_iTest
