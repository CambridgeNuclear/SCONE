program eigenvalue

  use numPrecision
  use genericProcedures,          only : printStart
  use IOdictionary_class,         only : IOdictionary
  use physicsPackage_inter,       only : physicsPackage
  use physicsPackageFactory_func, only : new_physicsPackage

  implicit none
  type(IOdictionary)                :: input
  class(physicsPackage),allocatable :: core

  call printStart()

  call input % initFrom('./InputFiles/FirstInput.c')

  allocate( core, source = new_physicsPackage(input))

  ! Read data

  call core % init(input)

  call core % run()


end program eigenvalue
