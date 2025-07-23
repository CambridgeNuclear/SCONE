program scone

  use numPrecision
  use display_func,               only : printStart, statusMsg
  use openmp_func,                only : ompSetNumThreads
  use mpi_func,                   only : mpiInit, mpiFinalise
  use commandLineUI,              only : getInputFile, clOptionIsPresent, addClOption, getFromCL
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
  integer(shortInt)                 :: cores

  ! Add command line options here
  call addClOption('--plot', 0, ['int'],&
          'Executes geometry plotting specified by a viz dict in the input file')
#ifdef _OPENMP
  call addClOption('--omp', 1, ['int'], &
          'Number of OpenMP threads in a parallel calculation')
#endif

  ! Get path to input file
  call getInputFile(inputPath)

  ! Initialize MPI
  call mpiInit()

  ! Set Number of threads
  if (clOptionIsPresent('--omp')) then
    call getFromCL(cores, '--omp', 1)
  else
    cores = 1
  end if
  call ompSetNumThreads(cores)

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

  call mpiFinalise()

  print *, 'Total calculation time: ', trim(secToChar(timerTime(timerIdx)))
  print *, 'Have a good day and enjoy your results analysis!'

end program scone
