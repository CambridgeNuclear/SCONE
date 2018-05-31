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
  type(eigenPhysicsPackage :: core

  ! Read data
  call matData   % initFrom('./materialInput')
  call geomData  % initFrom('./lattice.txt')
  call transData % initFrom('./transOp.txt')
  call collData  % initFrom('./collOp.txt')
  call inactiveTally % initFrom('./iaTally.txt')
  call activeTally   % initFrom('./aTally.txt')

  call core % init(matData, geomData, collData, transData, inactiveTally, activeTally)

  call core % run()


end program eigenvalue
