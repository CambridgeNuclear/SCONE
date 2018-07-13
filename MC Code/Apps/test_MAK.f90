program test

  use numPrecision
  use RNG_class
  use genericProcedures
  !use hashFunctions_func, only : FNV_1
  !use vector_class, only : vector
  use surface_inter, only : surfaceShelf, surface
  use sphere_class,  only : sphere

  implicit none



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


