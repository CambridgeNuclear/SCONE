program test
  
  use numPrecision
  use dictionary_class
  use IOdictionary_class
  use nuclearDataRegistry_mod

  implicit none

  type(IOdictionary) :: dict

  call dict % initFrom('./InputFiles/FullLib.c')

  call build_NuclearData(dict % getDictPtr('nuclearData'))

  call dict % kill()
  !call kill_NuclearData()

end program test


