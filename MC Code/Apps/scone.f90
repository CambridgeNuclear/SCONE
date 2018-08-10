program scone

  use numPrecision
  use genericProcedures,          only : printStart
  use commandLineUI,              only : getInputFile
  use IOdictionary_class,         only : IOdictionary
  use physicsPackage_inter,       only : physicsPackage
  use physicsPackageFactory_func, only : new_physicsPackage

  implicit none
  type(IOdictionary)                :: input
  class(physicsPackage),allocatable :: core
  character(:),allocatable          :: inputPath

  ! Add command line options here

  ! Get path to input file
  call getInputFile(inputPath)

  call printStart()

  call input % initFrom(inputPath)

  allocate( core, source = new_physicsPackage(input))

  call core % run()

end program scone
