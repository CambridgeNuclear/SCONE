!!
!! Module to help with simplify the process of writing input files, specifically aimed at IMC but
!! may be useful for other applications if needed
!!
module discretiseGeom_class

  use numPrecision
  use universalVariables
  use genericProcedures,              only : fatalError, numToChar
  use dictionary_class,               only : dictionary
  use dictParser_func,                only : fileToDict

  ! Geometry
  use geometry_inter,                 only : geometry
  use geometryReg_mod,                only : gr_geomPtr  => geomPtr, &
                                             gr_addGeom => addGeom, &
                                             gr_geomIdx  => geomIdx, &
                                             gr_kill     => kill
  ! Nuclear Data
  use materialMenu_mod,               only : mm_matTemp        => matTemp ,&
                                             mm_matFile        => matFile
  use nuclearDataReg_mod,             only : ndReg_init        => init ,&
                                             ndReg_activate    => activate ,&
                                             ndReg_kill        => kill, &
                                             ndReg_get         => get
  use nuclearDatabase_inter,          only : nuclearDatabase

  ! Useful variables
  class(geometry), pointer                     :: inputGeom
  class(nuclearDatabase), pointer              :: inputNucData
  integer(shortInt), dimension(:), allocatable :: sizeN

contains

  !!
  !! Initialises geometry and nuclear database given in input file, and uses this as a baseline to
  !! write new geometry and nuclear data dictionaries by splitting the input geometry into a
  !! uniform grid. Makes new text files 'newGeom.txt' and 'newData.txt'.
  !!
  !! New geometry is a lattice universe, with each lattice cell containing the same material at its
  !! centre in the initial input geometry, but with each cell initialising a new instance of this
  !! material object, allowing for spatial temperature variation.
  !!
  !! For N lattice cells, will write to files N instances of:
  !!   In newGeom.txt =>    pj { id j; type pinUniverse; radii (0); fills (mj); }
  !!   In newData.txt =>    mj { temp ...; composition {} xsFile ./.../...; volume ...; }
  !!                                          (filled in with correct data for underlying material)
  !!
  !! In future, might be able to somehow modify the way that geometry and materials are built such
  !! that separate instances of each mat are not needed
  !!
  !! Sample input:
  !!   dicretise { dimensions (10 1 1); }
  !!
  !! Args:
  !!   dict [in]     -> Input file dictionary
  !!   newGeom [out] -> dictionary containing new geometry input
  !!   newData [out] -> dictionary containing new nuclear database input
  !!
  !! Errors:
  !!   invalid input dimensions
  !!
  subroutine discretise(dict, newGeom, newData)
    class(dictionary), intent(in)                :: dict
    type(dictionary), intent(out)                :: newGeom
    type(dictionary), intent(out)                :: newData
    class(dictionary), pointer                   :: tempDict
    character(nameLen)                           :: dataName, geomName
    integer(shortInt)                            :: inputGeomIdx
    character(100), parameter :: Here = 'discretise (discretiseGeom_class.f90)'

    ! Get discretisation settings
    tempDict => dict % getDictPtr('discretise')
    call tempDict % get(sizeN, 'dimensions')
    if (any(sizeN <= 0)) call fatalError(Here, 'Invalid dimensions given')

    ! Build Nuclear Data using input
    call ndReg_init(dict % getDictPtr("nuclearData"))

    ! Build geometry using input
    tempDict => dict % getDictPtr('geometry')
    geomName = 'inputGeom'
    call gr_addGeom(geomName, tempDict)
    inputGeomIdx = gr_geomIdx(geomName)
    inputGeom    => gr_geomPtr(inputGeomIdx)

    ! Activate input Nuclear Data
    dataName = 'mg'
    call ndReg_activate(P_PHOTON_MG, dataName, inputGeom % activeMats())
    inputNucData => ndReg_get(P_PHOTON_MG)

    ! Create files for materials and geometry
    open(unit = 21, file = 'newGeom.txt')
    open(unit = 22, file = 'newData.txt')

    ! Write required information to files
    call writeToFiles(tempDict)

    close(21)
    close(22)

    ! Kill input geom and nuclear database
    call inputGeom    % kill()
    call inputNucData % kill()
    call gr_kill()
    call ndReg_kill()

    ! Create geometry and nuclear database from new files
    call fileToDict(newData, './newData.txt')
    call fileToDict(newGeom, './newGeom.txt')

  end subroutine discretise

  !!
  !! Automatically write required information to files in the following format:
  !!
  !! newGeom.txt:
  !!   type geometryStd;
  !!   boundary (...);
  !!   graph {type shrunk;}
  !!   surfaces { outer {id 1; type box; origin (...); halfwidth (...);}}
  !!   cells {}
  !!   universes {
  !!   root {id rootID; type rootUniverse; border 1; fill u<latID>;}
  !!   lattice {id latID; type latUniverse; origin (0 0 0); pitch (...); shape (nx ny nz); padMat void;
  !!   map (
  !!   1
  !!   2
  !!   .
  !!   .
  !!   .
  !!   nx*ny*nz); }
  !!   p1 {id 1; type pinUniverse; radii (0); fills (m1);}
  !!   p2 {id 2; type pinUniverse; radii (0); fills (m2);}
  !!   .
  !!   .
  !!   .
  !!   pN {id N; type pinUniverse; radii (0); fills (mN);} }
  !!
  !! newData.txt:
  !!   handles {mg {type baseMgIMCDatabase;}}
  !!   materials {
  !!   m1 {temp ...; composition {} xsFile ./...; volume ...;}
  !!   m2 {temp ...; composition {} xsFile ./...; volume ...;}
  !!   .
  !!   .
  !!   .
  !!   mN {temp ...; compositino {} xsFile ./...; volume ...;} }
  !!
  !!
  !! Args:
  !!   dict [in] -> geometry input dictionary
  !!
  subroutine writeToFiles(dict)
    class(dictionary), intent(in), pointer       :: dict
    character(nameLen)                           :: geomType, graphType
    integer(shortInt), dimension(:), allocatable :: boundary
    class(dictionary), pointer                   :: tempDict
    real(defReal)                                :: volume
    real(defReal), dimension(3)                  :: origin, halfwidth, pitch, corner, r
    real(defReal), dimension(6)                  :: bounds
    integer(shortInt)                            :: i, N, rootID, latID, matIdx, uniqueID
    integer(shortInt), dimension(3)              :: ijk

    ! Get bounds, pitch and lower corner
    bounds = inputGeom % bounds()
    do i = 1, 3
      pitch(i)  = (bounds(i+3)-bounds(i)) / sizeN(i)
      corner(i) = bounds(i)
    end do

    ! Calculate volume of each cell and number of cells
    volume = pitch(1) * pitch(2) * pitch(3)
    N = sizeN(1) * sizeN(2) * sizeN(3)

    ! Get geometry settings from input file
    call dict % get(geomType, 'type')
    call dict % get(boundary, 'boundary')
    tempDict => dict % getDictPtr('graph')
    call tempDict % get(graphType, 'type')

    ! Write to new file
    write (21, '(8A)') 'type '//trim(geomType)//';'//new_line('A')//&
                       &'boundary ('//numToChar(boundary)//');'//new_line('A')//&
                       &'graph {type '//trim(graphType)//';}'

    ! Construct outer boundary surface based on input geometry bounds
    ! Will likely cause problems if input outer surface is not a box
    do i=1, 3
      origin(i)    = (bounds(i+3) + bounds(i)) / 2
      halfwidth(i) = (bounds(i+3) - bounds(i)) / 2
    end do
    write (21, '(8A)') 'surfaces { outer {id 1; type box; origin ('//numToChar(origin)//'); &
                       &halfwidth ('//numToChar(halfwidth)//');}}'

    ! Empty cell dictionary and construct root universe
    rootID = N + 1
    latID  = N + 2
    write (21, '(8A)') 'cells {}'//new_line('A')//&
                       &'universes {'//new_line('A')//&
                       &'root {id '//numToChar(rootID)//'; type rootUniverse; border 1; &
                       &fill u<'//numToChar(latID)//'>;}'

    ! Write lattice input
    write(21, '(8A)') 'lattice {id '//numToChar(latID)//'; type latUniverse; origin (0 0 0); pitch ('//&
                       &numToChar(pitch)//'); shape ('//numToChar(sizeN)//'); padMat void;'&
                       &//new_line('A')//'map ('
    do i = 1, N
      write(21, '(8A)') numToChar(i)
    end do

    write(21, '(8A)') '); }'//new_line('A')


    ! Write material settings - TODO obtain this from input file automatically?
    write (22, '(8A)') 'handles {mg {type baseMgIMCDatabase;}}'//new_line('A')//'materials {'


    ! Write pin universes and materials
    do i = 1, N

      ! Get location in centre of lattice cell
      ijk = get_ijk(i, sizeN)
      r = corner + pitch * (ijk - HALF) 

      ! Get matIdx from input geometry
      call inputGeom % whatIsAt(matIdx, uniqueId, r)

      ! Pin universe for void and outside mat
      if (matIdx == VOID_MAT) then
        write(21, '(8A)') 'p'//numToChar(i)//' {id '//numToChar(i)//'; type pinUniverse; &
                          &radii (0); fills (void);}'

      else if (matIdx == OUTSIDE_MAT) then
        write(21, '(8A)') 'p'//numToChar(i)//' {id '//numToChar(i)//'; type pinUniverse; &
                          &radii (0); fills (outside);}'

      ! User-defined materials
      else
        ! Pin universe
        write(21, '(8A)') 'p'//numToChar(i)//' {id '//numToChar(i)//'; type pinUniverse; &
                          &radii (0); fills (m'//numToChar(i)//');}'

        ! Material
        write(22, '(8A)') 'm'//numToChar(i)//' {temp '//numToChar(mm_matTemp(matIdx))//'; &
                          &composition {} xsFile '//trim(mm_matFile(matIdx))//'; volume '//numToChar(volume)//';}'

      end if
    end do

    ! Write closing brackets for dictionaries
    write (21, '(8A)') '}'
    write (22, '(8A)') '}'

  end subroutine writeToFiles

  !!
  !! Generate ijk from localID and shape
  !!
  !! Args:
  !!   localID [in] -> Local id of the cell between 1 and product(sizeN)
  !!   sizeN [in]   -> Number of cells in each cardinal direction x, y & z
  !!
  !! Result:
  !!   Array ijk which has integer position in each cardinal direction
  !!
  pure function get_ijk(localID, sizeN) result(ijk)
    integer(shortInt), intent(in)               :: localID
    integer(shortInt), dimension(3), intent(in) :: sizeN
    integer(shortInt), dimension(3)             :: ijk
    integer(shortInt)                           :: temp, base

    temp = localID - 1

    base = temp / sizeN(1)
    ijk(1) = temp - sizeN(1) * base + 1

    temp = base
    base = temp / sizeN(2)
    ijk(2) = temp - sizeN(2) * base + 1

    ijk(3) = base + 1

  end function get_ijk


end module discretiseGeom_class
