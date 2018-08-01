program geomPrint

  use numPrecision
  use RNG_class
  use genericProcedures
  use universalVariables

  use csg_class,               only : csg
  use basicCellCSG_class,      only : basicCellCSG
  use IOdictionary_class,      only : IOdictionary
  use datalessMaterials_class, only : datalessMaterials


  implicit none
  type(basicCellCSG)                           :: geom
  type(datalessMaterials)                      :: nucData
  type(IOdictionary)                           :: geomData
  integer(shortInt)                            :: i, file
  integer(shortInt),dimension(:,:),allocatable :: colorMat

  call nucData % init(['uo2  ','uo211','uo222','uo233','water'])

  call geomData  % initFrom('./InputFiles/pinCell5.txt')

  call geom % init(geomData, nucData)

  !*** SLICE PLOT
  allocate(colorMat(1500,1500))
  call geom % slicePlot(colorMat,[ZERO, ZERO, ZERO],'z','material',[2.1_8,2.1_8])

  ! Print matrix to MATLAB Pcolor

  file = 7
  open(file, file='picture.m', action = 'write')

  write(file,*) "a = ["
  do i=1,size(colorMat,2)-1
     write(file,*) colorMat(:,i), ';'

  end do
  write(file,*) colorMat(:,i), '];'
  write(file,*) "figure"
  write(file,*) "h = pcolor(a)"
  write(file,*) "set(h,'EdgeColor','none')"
  write(file,*) "colorbar"

end program geomPrint


