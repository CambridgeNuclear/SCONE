program test_init_geom
  use numPrecision
  use genericProcedures
  use universalVariables
  use surface_inter
  use cell_class
  use universe_class
  use lattice_class
  use geometry_class
  use dictionary_class
  use IOdictionary_class
  use datalessMaterials_class, only : datalessMaterials
  use particle_class
  use rng_class
  use transportOperator_inter
  use transportOperatorDT_class
  use transportOperatorST_class
  use mesh_class
  use outputMesh_class

  implicit none

  type(geometry),pointer :: geom
  integer(shortInt) :: max_nest
  class(surface_ptr), dimension(:), allocatable :: surfaces
  type(IOdictionary) :: geomDict
  type(datalessMaterials) :: materials
  real(defReal) :: pitch = 1.0
  real(defReal), dimension(3) :: orgGeom
  integer(shortInt) :: pixels = 1000
  integer(shortInt) :: i,j
  integer(shortInt), dimension(:,:), allocatable :: colourMatrix
  character(100) :: path = './InputFiles/pinCell.txt'
  type(particle) :: p
  type(rng), target :: rand
  type(cell_ptr) :: c
  type(transportOperatorDT) :: DTOperator
  type(transportOperatorST) :: STOperator
  type(mesh) :: m
  real(defReal) :: angle
  type(outputMesh) :: vtk
  integer(shortInt), dimension(400,400,1) :: matMatrix
  integer(longInt), dimension(400,400,1) :: d

  call geomDict % initFrom(path)

  call materials % init(['uo2  ','water','big  '])
  print *, materials % materials

  !call initGeometryFromDict(geom, surfaces, max_nest, geomDict, materials)
  allocate(geom)
  call geom % init(geomDict, materials)

  print *,'Plotting geometry'
  allocate(colourMatrix(pixels,pixels))
  !colourMatrix = geom % slicePlot([pixels,pixels],orgGeom,[4.0*pitch,4.0*pitch,2.0*pitch],3)

  do i=1,pixels
    do j=1,pixels
      !print *,colourMatrix(i,j)
    end do
  end do

  print *,'Preparing transport'
  call p % build([0.2_8,0.0_8,ZERO],[sin(PI*0.9),cos(PI*0.9),ZERO],1,1._8)
  c = geom % whichCell(p % coords)
  call rand % init(500_8)

  ! New initialisation procedures require nuclear data
  !call DTOperator % init(geom)
  !call STOperator % init(geom)
  ! Circumvent the problem by manualy assigning pointer
  DTOperator % geom => geom
  STOperator % geom => geom

  call geom % constructBoundingBox([ZERO,ZERO,ZERO],[2._defReal,2._defReal,HALF])

  print *,'Calculating volumes'
 ! call geom % calculateVolumes(1000000,1_longInt)
  print *,'Finished volume calculation'
  print *,'Performing voxel calculation'
  call vtk % init([-2.0_defReal,-2.0_defReal,-HALF],[0.011_defReal, 0.011_defReal, ONE] ,[400,400,1])
  matMatrix = geom % voxelPlot(vtk % nVox, vtk % corner, vtk % width)
  call vtk % addData(real(matMatrix,defReal), 'materials')
  print *,'Finished voxel plot'
  call vtk % output('testVTK')
  stop


  ! Initialise mesh to score points
  call m % init([0.011_8,0.011_8,1._8], [-2.2_defReal,-2.2_defReal,-HALF], [400,400,1], 1, .FALSE.,'m1')

  print *, 'Begin walk'
  do i=1,1000000
    !f (modulo(i,1000) == 0) print *, i
    angle=2*PI*rand % get()
    call p % build([0.0_8,0.0_8,ZERO],[cos(3.0*PI/4),sin(3.0*PI/4),ZERO],1,1._8)
    p % pRNG => rand
    c = geom % whichCell(p % coords)
    !print *,c%name()
    do while (.not. p % isDead)
      call DTOperator % transport(p)
      call m % score(p % rGlobal(), p % dirGlobal())
      if (rand % get() > 0.95) then
        p % isDead = .TRUE.
      else
        angle=2*PI*rand % get()
        call p % point([cos(angle),sin(angle),ZERO])
      end if
    end do
  end do
  print *,'walk done'

  d = m % density
  call vtk % addData(real(d,8), 'particleDensity')
  call vtk % output('testVTK')

  !do i=1,200
  !  do j=1,200
  !    print *, m % density(i,j,1)
  !  end do
  !end do

end program test_init_geom
