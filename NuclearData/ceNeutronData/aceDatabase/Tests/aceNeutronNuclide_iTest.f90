module aceNeutronNuclide_iTest

  use numPrecision
  use aceCard_class,           only : aceCard
  use ceNeutronDatabase_inter, only : ceNeutronDatabase
  use aceNeutronNuclide_class, only : aceNeutronNuclide
  use pfUnit_mod

  implicit none


contains

  !!
  !! Load and verify functions of aceNeutronNuclide initialised to O-16
  !!
@Test
  subroutine testACEnuclideO16()
    type(aceNeutronNuclide), target   :: nuc
    type(aceCard)                     :: ACE
    class(ceNeutronDatabase), pointer :: database => null()

    ! Build ACE library
    call ACE % readFromFile('./IntegrationTestFiles/8016JEF311.ace', 1)

    ! Build nuclide
    call nuc % init(ACE, 1, database)

  end subroutine testACEnuclideO16


  !!
  !! Load and verify functions of aceNeutronNuclide initialised to U-233
  !!
@Test
  subroutine testACEnuclideU233()
    type(aceNeutronNuclide), target   :: nuc
    type(aceCard)                     :: ACE
    class(ceNeutronDatabase), pointer :: database => null()

    ! Build ACE library
    call ACE % readFromFile('./IntegrationTestFiles/92233JEF311.ace', 1)

    ! Build nuclide
    call nuc % init(ACE, 1, database)

  end subroutine testACEnuclideU233


end module aceNeutronNuclide_iTest
