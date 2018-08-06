program dynamic

  use numPrecision

  use IOdictionary_class,        only : IOdictionary
  use dynamPhysicsPackage_class, only : dynamPhysicsPackage

  implicit none

  type(IOdictionary)        :: matData
  type(IOdictionary)        :: geomData
  type(IOdictionary)        :: transData
  type(IOdictionary)        :: collData
  type(IOdictionary)        :: timeTally
  !type(IOdictionary)        :: inactiveTally
  type(dynamPhysicsPackage) :: core

  ! Read data
  call matData   % initFrom('./InputFiles/materialInput')
  call geomData  % initFrom('./InputFiles/pinCell.txt')
  call transData % initFrom('./InputFiles/transOp.txt')
  call collData  % initFrom('./InputFiles/collOp.txt')
  !call inactiveTally % initFrom('./InputFiles/iaTally.txt')
  call timeTally   % initFrom('./InputFiles/tTally.txt')

  call core % init(matData, geomData, collData, transData, timeTally)

  call core % run()


end program dynamic
