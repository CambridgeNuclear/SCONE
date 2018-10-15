program test
  use numPrecision

  use nuclearData_inter, only : nuclearData
  use IOdictionary_class, only : IOdictionary
  use nuclearDataRegistry_mod


  implicit none
  class(nuclearData),pointer :: myPtr => null()
  type(IOdictionary)  :: inDict
  class(dictionary),pointer :: dict
  character(nameLen),dimension(:),allocatable :: matNames
  character(nameLen)                          :: type
  integer(shortInt) :: i

  call inDict % initFrom('./InputFiles/FirstInput.c')

  dict => inDict % getDictPtr('materials')
  call dict % keysDict(matNames)
  type = 'transMG'


  myPtr => new_nuclearData_ptr(dict,type,matNames)

  do i=1,size(matNames)
    print *, matNames(i), 'INSIDE: ' , myPtr % getName(i)
  end do




end program test


