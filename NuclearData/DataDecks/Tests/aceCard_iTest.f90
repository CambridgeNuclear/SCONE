module aceCard_iTest

  use numPrecision
  use aceCard_class, only : aceCard
  use pFUnit_mod

  implicit none


contains


  !!
  !! Test Reloading aceCard
  !!
  !! Test Correcness of fission-related data after reloading
  !!
@Test
  subroutine testAceReloading_fissionData()
    type(aceCard) :: ACE
    character(*),parameter :: path1 = './IntegrationTestFiles/91231JEF311.ace'
    character(*),parameter :: path2 = './IntegrationTestFiles/91232JEF311.ace'
    real(defReal), parameter :: TOL = 1.0E-6_defReal

    ! Load 1st Card with Delayed Fission Neutrons
    call ACE % readFromFile(path1, 1)

    ! Test HEADER contents
    ! Mass
    @assertEqual(229.05_defReal,     ACE % AW, TOL* ACE % AW )
    ! Temperarture
    @assertEqual(2.5852E-08_defReal, ACE % TZ, TOL * ACE % TZ)
    @assertEqual('91231.03c', ACE % ZAID)

    ! Test Fission Data
    @assertEqual(8, ACE % precursorGroups())
    @assertTrue(ACE % isFissile())
    @assertTrue(ACE % hasNuPrompt())
    @assertTrue(ACE % hasNuDelayed())
    @assertTrue(ACE % hasNuTotal())

    ! Load 2nd Card without Delayed Fission Neutrons
    call ACE % readFromFile(path2, 1)

    ! Test HEADER Contents
    ! Mass
    @assertEqual(230.045_defReal,    ACE % AW, TOL* ACE % AW )
    ! Temperature
    @assertEqual(2.5852E-08_defReal, ACE % TZ, TOL * ACE % TZ)
    @assertEqual('91232.03c', ACE % ZAID)

    ! Test Fission Data
    @assertEqual(0, ACE % precursorGroups())
    @assertTrue(ACE % isFissile())
    @assertTrue(ACE % hasNuPrompt())
    @assertFalse(ACE % hasNuDelayed())
    @assertTrue(ACE % hasNuTotal())

  end subroutine testAceReloading_fissionData



end module aceCard_iTest
