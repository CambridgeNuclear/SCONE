program scone

  use numPrecision
  use genericProcedures,          only : printStart
  use commandLineUI,              only : getInputFile, clOptionIsPresent, addClOption
  use dictionary_class,           only : dictionary
  use dictParser_func,            only : fileToDict
  use physicsPackage_inter,       only : physicsPackage
  use physicsPackageFactory_func, only : new_physicsPackage
  use vizPhysicsPackage_class,    only : vizPhysicsPackage
  use timer_mod                 , only : registerTimer, timerStart, timerStop, timerTime, secToChar

  implicit none
  type(dictionary)                  :: input
  class(physicsPackage),allocatable :: core
  character(:),allocatable          :: inputPath
  integer(shortInt)                 :: timerIdx

  ! Add command line options here
  call addClOption('--plot',0,['int'],&
          'Executes geometry plotting specified by a viz dict in the input file')

  ! Get path to input file
  call getInputFile(inputPath)

  ! Register timer
  timerIdx = registerTimer('Main Timer')

  call printStart()

  call fileToDict(input, inputPath)

  call timerStart(timerIdx)

  if (clOptionIsPresent('--plot')) then
    allocate(vizPhysicsPackage :: core)
    call core % init(input)
  else
    allocate( core, source = new_physicsPackage(input))
  endif

  call core % run()

  call timerStop(timerIdx)
  print *, 'Total calculation time: ', trim(secToChar(timerTime(timerIdx)))
  print *, 'Have a good day and enjoy your result analysis!'
end program scone
