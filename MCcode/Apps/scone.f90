program scone

  use numPrecision
  use genericProcedures,          only : printStart
  use commandLineUI,              only : getInputFile
  use IOdictionary_class,         only : IOdictionary
  use physicsPackage_inter,       only : physicsPackage
  use physicsPackageFactory_func, only : new_physicsPackage
  use timer_mod                 , only : registerTimer, timerStart, timerStop, timerTime, secToChar

  implicit none
  type(IOdictionary)                :: input
  class(physicsPackage),allocatable :: core
  character(:),allocatable          :: inputPath
  integer(shortInt)                 :: timerIdx

  ! Add command line options here

  ! Get path to input file
  call getInputFile(inputPath)

  ! Register timer
  timerIdx = registerTimer('Main Timer')

  call printStart()

  call input % initFrom(inputPath)

  call timerStart(timerIdx)

  allocate( core, source = new_physicsPackage(input))

  call core % run()

  call timerStop(timerIdx)
  print *, 'Total calculation time: ', trim(secToChar(timerTime(timerIdx)))
  print *, 'Have a good day and enjoy your result analysis!'
end program scone
