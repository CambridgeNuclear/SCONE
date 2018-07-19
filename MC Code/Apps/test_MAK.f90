program test

  use numPrecision
  use RNG_class
  use genericProcedures
  use universalVariables
  !use hashFunctions_func, only : FNV_1
  use vector_class, only : vector
  use surface_inter, only : surfaceShelf, surface
  use sphere_class,  only : sphere
  use cell_class,    only : cell, cellShelf

  implicit none
  type(vector)               :: r, u
  real(defReal)              :: dist
  class(surface),allocatable :: s1
  class(surface),allocatable :: s2
  class(surface),allocatable :: s3
  type(surfaceShelf) :: sShelf
  type(cell)         :: c1
  type(cell)         :: c2
  type(cellShelf)    :: cShelf
  integer(shortInt)  :: i,j
  integer(shortInt),dimension(:),allocatable :: intArr


  allocate(sphere :: s1)
  allocate(sphere :: s2)
  allocate(sphere :: s3)


  select type(s1)
    type is (sphere)
      call s1 % init([ZERO, ZERO, ZERO], 1.0_8, 1)
  end select

  select type(s2)
    type is (sphere)
      call s2 % init([ZERO, ZERO, ZERO], 2.0_8, 2)
  end select

  select type(s3)
  type is (sphere)
      call s3 % init([ZERO, ZERO, ZERO], 3.0_8, 3)
  end select


  call sShelf % init(1)
  call sShelf % addUnique(s1,i)
  call sShelf % addUnique(s2,i)
  call sShelf % addUnique(s3,i)

  call c1 % init([-1],1,sShelf)
  call c2 % init([1, -3, -1],2,sShelf)


  call cShelf % init(1,1)
  call cShelf % addUnique(c1,i)
  call cShelf % addUnique(c2,i)

  r = [1.0_8 + 0.9_8 *surface_tol,0.0_8,0.0_8]
  u = [1.0,0.0,0.0]

  print *, cShelf % shelf(2) % isInside(r,u,sShelf)
  call cShelf % shelf(2) % distance(dist,i,r,u,sShelf)
  print *, dist

  call c2 % init([1, -3],7,sShelf)
  call cShelf % getOrAdd(c2,i,j)
  print *, i,j

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


