program test
  
  use numPrecision
  use nuclearData_inter,          only : nuclearData
  use transportNuclearData_inter, only : transportNuclearData

  use testTransportNuclearData_class,      only : testTransportNuclearData
  use byNucMT_class,              only : byNucMT


  implicit none
  class(nuclearData),pointer           :: xsDat   => null()
  class(transportNuclearData),pointer  :: transND => null()

  type(testTransportNuclearData),pointer :: target1
  type(byNucMT),pointer         :: target2


  allocate(target1)
  allocate(target2)

 ! call ptr(target1, xsDat)




contains

  subroutine ptr(a1, a2)
    class(nuclearData),pointer :: a1
    class(nuclearData),pointer :: a2

    print *, extends_type_of(a2, a1)
    a1 => a2

  end subroutine


end program test




