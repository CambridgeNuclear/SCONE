program test
  use numPrecision

  use nuclearData_inter, only : nuclearData
  use IOdictionary_class, only : IOdictionary
  use dictionary_class,   only : dictionary
  use nuclearDataRegistry_mod, only : new_nuclearData_ptr, nuclearData_buildMaterials, nuclearData_kill

  implicit none
  type(IOdictionary)  :: inDict
  class(dictionary),pointer :: dict
  character(nameLen),dimension(:),allocatable :: matNames


  call inDict % initFrom('./InputFiles/FirstInput.c')
  dict => inDict % getDictPtr('nuclearData')
  call nuclearData_buildMaterials(dict)
  call nuclearData_kill()
  call inDict % kill()

end program test


