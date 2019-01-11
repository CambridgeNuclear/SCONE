program test
  
  use numPrecision
  use tallyMap_inter,       only : tallyMap
  use tallyMapSlot_class,   only : tallyMapSlot
  use tallyMapFactory_func, only : new_tallyMap
  use dictionary_class,     only : dictionary
  
  implicit none
  type(tallyMapSlot) :: slot
  type(dictionary)   :: tempDict
  integer(shortInt)  :: i

  do i=1,100

  ! Build map lin
  call tempDict % init(5)
  call tempDict % store('type','energyMap')
  call tempDict % store('grid','lin')
  call tempDict % store('min', 0.01_defReal)
  call tempDict % store('max', 10.0_defReal)
  call tempDict % store('N', 20)

  call slot % init(tempDict)

  call tempDict % kill()

  end do

  call slot % kill()

end program test


