module elasticScattering_iTest

  use pFUnit_mod
  implicit none

contains

@Test
  subroutine myTest()

    @assertEqual(1,2)

  end subroutine myTest


    
end module elasticScattering_iTest
