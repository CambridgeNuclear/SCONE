program eigenvalue

  use numPrecision

  use IOdictionary_class,        only : IOdictionary
  use eigenPhysicsPackage_class, only : eigenPhysicsPackage

  implicit none

  type(IOdictionary)       :: matData
  type(IOdictionary)       :: geomData
  type(IOdictionary)       :: transData
  type(IOdictionary)       :: collData
  type(IOdictionary)       :: activeTally
  type(IOdictionary)       :: inactiveTally
  type(eigenPhysicsPackage) :: core

  ! Read data
  call matData   % initFrom('./InputFiles/materialInput')
  call geomData  % initFrom('./InputFiles/pinCell2.txt')
  call transData % initFrom('./InputFiles/transOp.txt')
  call collData  % initFrom('./InputFiles/collOp.txt')
  call inactiveTally % initFrom('./InputFiles/iaTally.txt')
  call activeTally   % initFrom('./InputFiles/aTally.txt')

  call core % init(matData, geomData, collData, transData, inactiveTally, activeTally)

  call core % run()


end program eigenvalue
