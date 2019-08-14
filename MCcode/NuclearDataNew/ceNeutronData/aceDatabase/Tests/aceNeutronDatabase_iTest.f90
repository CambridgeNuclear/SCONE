module aceNeutronDatabase_iTest

  use numPrecision
  use dictionary_class,         only : dictionary
  use IOdictionary_class,       only : IOdictionary
  use aceNeutronDatabase_class, only : aceNeutronDatabase
  use nuclearDatabase_inter,    only : nuclearDatabase
  use materialMenu_mod,         only : mm_init => init, mm_kill => kill
  use pFUnit_mod

  implicit none

  ! Material definitions
  character(*),parameter :: MAT_INPUT_STR = " & 
  water { temp 273;                           &
         composition {                        &
         1001.03 5.028E-02;                   &
         8016.03 2.505E-02;                   &
         }                                    &
       }                                      &
  uo2  { temp 1;                              &
          composition {                       &
          92233.03 2.286E-02;                 &
          8016.03  4.572E-02;                 &
          }                                   &
        }"

  ! CE Neutron Database specification
  character(*),parameter :: ACE_INPUT_STR = " & 
  aceLibrary ./IntegrationTestFiles/testLib;    " 

contains

  !!
  !! One big monster test to avoid expensive set up each test
  !!
@Test
  subroutine test_aceNeutronDatabase()
    type(aceNeutronDatabase), target :: data
    class(nuclearDatabase), pointer  :: ptr
    type(IOdictionary)               :: matDict
    type(IOdictionary)               :: dataDict

    ! Prepare dictionaries
    call matDict % initFromChar(MAT_INPUT_STR)
    call dataDict % initFromChar(ACE_INPUT_STR)

    ! Build material menu
    call mm_init(matDict)

    ! Initialise data
    ptr => data
    call data % init(dataDict, ptr, silent = .true.)


    ! Clean eeverything
    call mm_kill()

  end subroutine test_aceNeutronDatabase

end module aceNeutronDatabase_iTest
