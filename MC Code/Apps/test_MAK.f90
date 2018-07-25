program test

  use numPrecision
  use RNG_class
  use genericProcedures
  use universalVariables
  use hashFunctions_func, only : FNV_1, knuthHash
  use vector_class, only : vector
  use surface_inter, only : surfaceShelf, surface
  use sphere_class,  only : sphere
  use cell_class,    only : cell, cellShelf
  use maps_class,    only : intMap
  use cellUniverse_class, only : cellUniverse
  use coord_class,        only : coord

  use csg_class,          only : csg
  use basicCellCSG_class, only : basicCellCSG
  use IOdictionary_class, only : IOdictionary
  use datalessMaterials_class, only : datalessMaterials

  implicit none

  type(basicCellCSG) :: geom
  type(datalessMaterials) :: nucData
  type(IOdictionary) :: geomData
  integer(shortInt)      :: i, file
  integer(shortInt),dimension(:,:),allocatable :: colorMat

  call nucData % init(['uo2  ','water'])

  call geomData  % initFrom('./InputFiles/pinCell2.txt')

  call geom % init(geomData, nucData)

  allocate(colorMat(1500,1500))

  call geom % slicePlot(colorMat,[ZERO, ZERO, ZERO],'x','material',[5.0_8,5.0_8])

  ! Print matrix to MATLAB Pcolor

  file = 7
  open(file, file='picture.m', action = 'write')

  write(file,*) "a = ["
  do i=1,size(colorMat,2)-1
     write(file,*) modulo(colorMat(:,i),10), ';'

  end do
  write(file,*) modulo(colorMat(:,i),10), '];'
  write(file,*) "figure"
  write(file,*) "h = pcolor(a)"
  write(file,*) "set(h,'EdgeColor','none')"
  write(file,*) "colorbar"




!  type(vector)               :: r, u
!  real(defReal)              :: dist
!  class(surface),allocatable :: s1
!  class(surface),allocatable :: s2
!  class(surface),allocatable :: s3
!  type(surfaceShelf) :: sShelf
!  type(cell)         :: c1
!  type(cell)         :: c2
!  type(cell)         :: c3
!  type(cellShelf)    :: cShelf
!  type(cellUniverse) :: uni
!  type(coord)        :: coords
!  integer(shortInt)  :: i,j
!  integer(shortInt),dimension(:),allocatable :: intArr
!
!  allocate(sphere :: s1)
!  allocate(sphere :: s2)
!  allocate(sphere :: s3)
!
!
!  select type(s1)
!    type is (sphere)
!      call s1 % init([ZERO, ZERO, ZERO], 1.0_8, 1)
!  end select
!
!  select type(s2)
!    type is (sphere)
!      call s2 % init([ZERO, ZERO, ZERO], 2.0_8, 2)
!  end select
!
!  select type(s3)
!  type is (sphere)
!      call s3 % init([ZERO, ZERO, ZERO], 3.0_8, 3)
!  end select
!
!
!  call sShelf % init(1)
!  call sShelf % addUnique(s1,i)
!  call sShelf % addUnique(s2,i)
!  call sShelf % addUnique(s3,i)
!
!  call c1 % init([-1],2,sShelf)
!  call c2 % init([1, -3 ],1,sShelf)
!  call c3 % init([3],77,sShelf)
!
!  call cShelf % init(1,1)
!  call cShelf % add(c1,i)
!  call cShelf % add(c2,i)
!  call cShelf % add(c3,i)
!
!
!  call uni % init([ZERO, ZERO, ZERO],1,[2,77,1],cShelf)
!
!  coords % r   = [1.0_8- 0.9 * surface_tol, 0.0_8, 0.0_8]
!  coords % dir = [0.1, 0.6, 0.0]
!  coords % dir = coords % dir / norm2(coords % dir)
!
!  call uni % enter(coords, cShelf, sShelf)
!
!  call uni % distance(dist,i,coords,cShelf, sShelf)
!  print *, dist, i
!  coords % r = coords % r + coords % dir * dist
!
!  call uni % cross(coords, i, cShelf, sShelf)
!  call uni % distance(dist,i,coords,cShelf, sShelf)
!  print *, dist, i
!
!
!  print *, coords % r
!  print *, coords % dir
!  print *, coords % uniIdx, coords % uniRootID
!  print *, coords % localID, coords % cellIdx
!
!
!
!
!


!! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! intMap test
!! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!  type(intMap) :: myMap
!
!  call myMap % init(1)
!  call myMap % add(1,7)
!  call myMap % add(2,77)
!  call myMap % add(3,777)
!  call myMap % add(2,7777)
!  call myMap % add(4,77777)
!  call myMap % add(5,7)
!  call myMap % add(6,7)
!  !call myMap % add(7,7)
!  !call myMap % add(8,7)
!
!  print *, size(myMap % map), myMap % Nexp, myMap % N, myMap % Load
!  print *, myMap % map % key
!  print *, myMap % map % val
!
!  print *, myMap % get(1)
!  print *, myMap % get(2)
!  print *, myMap % get(1)
!  print *, myMap % get(3)
!  print *, myMap % get(4)
!  print *, myMap % get(99)


!! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! Cell & CellShelf tests
!! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!  type(vector)               :: r, u
!  real(defReal)              :: dist
!  class(surface),allocatable :: s1
!  class(surface),allocatable :: s2
!  class(surface),allocatable :: s3
!  type(surfaceShelf) :: sShelf
!  type(cell)         :: c1
!  type(cell)         :: c2
!  type(cellShelf)    :: cShelf
!  integer(shortInt)  :: i,j
!  integer(shortInt),dimension(:),allocatable :: intArr
!
!  allocate(sphere :: s1)
!  allocate(sphere :: s2)
!  allocate(sphere :: s3)
!
!
!  select type(s1)
!    type is (sphere)
!      call s1 % init([ZERO, ZERO, ZERO], 1.0_8, 1)
!  end select
!
!  select type(s2)
!    type is (sphere)
!      call s2 % init([ZERO, ZERO, ZERO], 2.0_8, 2)
!  end select
!
!  select type(s3)
!  type is (sphere)
!      call s3 % init([ZERO, ZERO, ZERO], 3.0_8, 3)
!  end select
!
!
!  call sShelf % init(1)
!  call sShelf % addUnique(s1,i)
!  call sShelf % addUnique(s2,i)
!  call sShelf % addUnique(s3,i)
!
!  call c1 % init([-1],1,sShelf)
!  call c2 % init([1, -3, -1],2,sShelf)
!
!
!  call cShelf % init(1,1)
!  call cShelf % addUnique(c1,i)
!  call cShelf % addUnique(c2,i)
!
!  r = [1.0_8 + 0.9_8 *surface_tol,0.0_8,0.0_8]
!  u = [1.0,0.0,0.0]
!
!  print *, cShelf % shelf(2) % isInside(r,u,sShelf)
!  call cShelf % shelf(2) % distance(dist,i,r,u,sShelf)
!  print *, dist
!
!  call c2 % init([1, -3],7,sShelf)
!  call cShelf % getOrAdd(c2,i,j)
!  print *, i,j

!! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!! SurfaceShelf tests
!! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!
!  type(surfaceShelf) :: surfBank
!  class(surface),allocatable       :: surf
!  character(:),allocatable :: def
!  integer(shortInt) :: sIdx, sId
!
!  allocate( sphere :: surf)
!  select type(surf)
!    type is (sphere)
!      call surf % init([ZERO, ZERO, ZERO], 5.0_8, 1)
!  end select
!
!
!  call surf % getDef(def)
!  print *, def
!  def = ''
!
!  call surfBank % init(1,100)
!
!  call surfBank % addUnique(surf,sIdx)
!  call surfBank % shelf(1) % getDef(def )
!  print *, def, sIdx
!
!  allocate( sphere :: surf)
!  select type(surf)
!    type is (sphere)
!      call surf % init([ZERO, ZERO, ZERO], 5.000001_8, 7)
!  end select
!
!
!  call surfBank % addUnique(surf,sIdx)
!  call surfBank % shelf(2) % getDef(def )
!  print *, def, sIdx
!
!
!
!  allocate( sphere :: surf)
!  select type(surf)
!    type is (sphere)
!      call surf % init([ZERO, ZERO, ZERO], 5.000003_8, 7)
!  end select
!
!  call surfBank % getOrAdd(surf,sId,sIdx)
!  print *, sId,sIdx
!
!  print *, size(surfBank % shelf)
!  call surfBank % trimSize()
!  print *, size(surfBank % shelf)
!
!  print *, surfBank % getIdx([1,2,3,4,5,6,7,8,9,10])

end program test


