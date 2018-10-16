program test
  use numPrecision

  use nuclearData_inter, only : nuclearData
  use IOdictionary_class, only : IOdictionary
  use dictionary_class,   only : dictionary
  use nuclearDataRegistry_mod, only : build_nuclearData, kill_nuclearData, getMatIdx, getHandlePtr

  implicit none
  type(IOdictionary)  :: inDict
  class(dictionary),pointer :: dict
  character(nameLen),dimension(:),allocatable :: matNames
  character(nameLen)                          :: name1, name2
  class(nuclearData), pointer                 :: nucDat

  call inDict % initFrom('./InputFiles/FirstInput.c')
  dict => inDict % getDictPtr('nuclearData')
  call build_nuclearData(dict)
  !call kill_nuclearData()

  name1 = 'uo2'
  name2 = 'water'

  print *, getMatIdx(name1), getMatIdx(name2)

  name1 = 'ce2'
  nucDat => getHandlePtr(name1)
  print *, nucDat % getName(1)

  call inDict % kill()

end program test


