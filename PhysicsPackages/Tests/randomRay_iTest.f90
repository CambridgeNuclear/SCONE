module randomRay_iTest

  use numPrecision
  use dictionary_class,  only : dictionary
  use dictParser_func,   only : charToDict
  use charMap_class,     only : charMap
  use dictParser_func,   only : fileToDict
  use randomRayPhysicsPackage_class, only : randomRayPhysicsPackage
  use pFUnit_mod

  implicit none


contains

  !!
  !! Random ray integration test.
  !! Homogeneous, reflected box. 
  !! Should have a precise eigenvalue, the ratio of nuSigmaF/SigmaA.
  !! Uncertainties will not be present as all points in phase space are identical.
  !!
@Test
  subroutine test_random_ray()
    type(randomRayPhysicsPackage) :: pp
    character(*), parameter       :: path = './IntegrationTestFiles/PhysicsPackages/SlabTRRM'
    type(dictionary)              :: dict
    real(defReal)                 :: keff
    real(defReal), parameter      :: TOL = 1.0E-5_defReal
    
    ! Load dictionary
    call fileToDict(dict, path)

    ! Initialise physics package and run
    call pp % init(dict, loud = .false.)
    call pp % run()

    ! Extract and verify keff
    keff = pp % keffScore(1)
    @assertEqual(keff, 2.25_defReal, TOL)

    ! Kill physics package
    call pp % kill()

  end subroutine test_random_ray

end module randomRay_iTest
