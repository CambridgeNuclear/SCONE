module renderPhysicsPackage_class

  use numPrecision
  use universalVariables
  use genericProcedures,              only : fatalError, rotateAroundAxis, crossProduct, numToChar
  use dictionary_class,               only : dictionary

  ! Physics package interface
  use physicsPackage_inter,           only : physicsPackage

  ! Geometry
  use geometry_inter,                 only : geometry
  use geometryReg_mod,                only : gr_geomPtr  => geomPtr, gr_geomIdx  => geomIdx
  use geometryFactory_func,           only : new_geometry

  ! Nuclear Data
  use materialMenu_mod,               only : mm_nMat => nMat, mm_matIdx => matIdx, mm_matName => matName
  use nuclearDataReg_mod,             only : ndReg_init        => init ,&
                                             ndReg_getMatNames => getMatNames
  use nuclearDatabase_inter,          only : nuclearDatabase

  ! Visualisation
  use imgBmp_func,                    only : imgBmp_toFile
  use visualiser_class,               only : visualiser, materialColour, brightnessScale, &
                                             getRayPlotInfo, getRayPlotTransformation

  implicit none
  private

  !!
  !! Physics Package for live rendering of ray plots
  !! As well as a geometry, requires only a visualiser dictionary
  !! containing a ray plot.
  !!
  !! The rendering allows modifying the ray plot live. Requires the presence
  !! of an app for viewing/updating images (eog works well).
  !!
  type, public,extends(physicsPackage) :: renderPhysicsPackage
    private
    ! Building blocks
    class(geometry), pointer :: geom => null()
    integer(shortInt)        :: geomIdx = 0
    type(visualiser)         :: viz
    type(dictionary)         :: dict
    character(nameLen)       :: displayApp = 'eog'

  contains
    procedure :: init
    procedure :: run
    procedure :: kill

  end type renderPhysicsPackage

contains

  !!
  !! Calls visualiser to generate visualisation
  !!
  subroutine run(self)
    class(renderPhysicsPackage), intent(inout) :: self
    real(defReal)                                  :: fov, ambient
    integer(shortInt), dimension(2)                :: res
    integer(shortInt), dimension(:), allocatable   :: mats, temp
    integer(shortInt), dimension(:,:), allocatable :: img, matIDs
    real(defReal), dimension(:,:), allocatable     :: lum
    real(defReal), dimension(3,3)                  :: M
    integer(shortInt)                              :: offset, ios, matIdx, pos, n, i, j
    character(nameLen)                             :: outputFile, outputLive
    real(defReal)                                  :: a, b, c, d, e, f
    real(defReal), dimension(3)                    :: dir, dirV, dirH, centre, camera,&
                                                      light, up, right
    real(defReal), dimension(6)                    :: bounds
    character(100)                                 :: line
    character(nameLen)                             :: cmd, matName, argStr
    character(100), parameter :: Here = 'run (renderPhysicsPackage_class.f90)'

    ! Extract ray info from the viz dictionary
    call getRayPlotInfo(self % dict, outputFile, centre, camera, light, up, M, fov, ambient, &
                        res, mats, offset, bounds)
    
    allocate(matIDs(res(1), res(2)))
    allocate(lum(res(1), res(2)))
      
    dirV = M(:,3)
    dirH = M(:,2)

    ! Overwrite output
    outputFile = 'output_temp.bmp'
    outputLive = 'output.bmp'
    
    print *, 'COMMANDS:'
    print *, '1D Pan camera: w/a/s/d <distance in cm>'
    print *, '2D Pan camera: pan <horizontal distance in cm> <vertical distance in cm>'
    print *, 'Rotate camera: rot <yaw in degrees> <pitch in degrees>'
    print *, 'Zoom camera: +/- <distance in cm>'
    print *, 'Set "up": up <3D direction>'
    print *, 'Set light: light <3D position>'
    print *, 'Set camera: camera <3D position>'
    print *, 'Set centre: centre <3D position>'
    print *, 'Set ambient light: amb <strength between 0 and 1>'
    print *, 'Make a material opaque: opaq <material name>'
    print *, 'Make a material transparent: transp <material name>'
    print *, 'Set rendered volume: bounds <-x -y -z +x +y +z>'
    print *, 'Set resolution in x and y directions: res <xPixels yPixels>'
    print *, 'Provide list of material names: mats'
    print *, 'Provide the positions of the camera, centre, and light: pos'
    print *, 'Quit: q'
      
    ! lum contains luminosity values, matIDs identifies which materials were hit
    call self % geom % rayPlot(lum, matIDs, camera, light, M, mats, fov, ambient, bounds)
    
    ! Translate to an image.
    ! Obtain material colours and scale by luminosity
    matIDs = materialColour(matIDs, offset)
    img = matIDs

    ! Scale image by brightness
    img = brightnessScale(img, lum)

    ! Print image
    call imgBmp_toFile(img, outputFile)
    call execute_command_line("mv "//outputFile//" "//outputLive)
    call execute_command_line(self % displayApp//" "//outputLive//" &")

    print *,'Input command:'
    do
    
      a = ZERO
      b = ZERO
      c = ZERO
      cmd = ''
      read(*,'(A)') line
      read(line, *) cmd

      pos = index(line, ' ')
      if (pos > 0) then
        argStr = adjustl(line(pos+1:))
      else
        argStr = ''
      end if
      
      select case (trim(cmd))
      case ("+")
        read(argStr, *, iostat=ios) a
        dir = centre - camera
        dir = dir / norm2(dir)
        camera = camera + a * dir
      case ("-")
        read(argStr, *, iostat=ios) a
        dir = centre - camera
        dir = dir / norm2(dir)
        camera = camera - a * dir
      case ("w")
        read(argStr, *, iostat=ios) a
        centre = centre + a * dirV
        camera = camera + a * dirV
      case ("s")
        read(argStr, *, iostat=ios) a
        centre = centre - a * dirV
        camera = camera - a * dirV
      case ("a")
        read(argStr, *, iostat=ios) a
        centre = centre - a * dirH
        camera = camera - a * dirH
      case ("d")
        read(argStr, *, iostat=ios) a
        centre = centre + a * dirH
        camera = camera + a * dirH
      case ("pan")
        read(argStr, *, iostat=ios) a, b
        centre = centre + a * dirH + b * dirV
        camera = camera + a * dirH + b * dirV
      case ("rot")
        read(argStr, *, iostat=ios) a, b
        
        dir = camera - centre
        
        ! Yaw
        dir = rotateAroundAxis(dir, up, a * PI / 180)
        
        ! Pitch
        right = crossProduct(dir, up)
        right = right / norm2(right)
        dir = rotateAroundAxis(dir, right, b * PI / 180)

        camera = centre + dir
        
        M = getRayPlotTransformation(camera, centre, up)
        dirV = M(:,3)
        dirH = M(:,2)

      case ("up")
        read(argStr, *, iostat=ios) a, b, c
        up = [a, b, c]
        if (all(up == [ZERO, ZERO, ZERO])) up = [ZERO, ZERO, ONE]
        up = up /norm2(up)
        M = getRayPlotTransformation(camera, centre, up)
        dirV = M(:,3)
        dirH = M(:,2)
      case ("camera")
        read(argStr, *, iostat=ios) a, b, c
        camera = [a, b, c]
        M = getRayPlotTransformation(camera, centre, up)
        dirV = M(:,3)
        dirH = M(:,2)
      case ("centre")
        read(argStr, *, iostat=ios) a, b, c
        centre = [a, b, c]
        M = getRayPlotTransformation(camera, centre, up)
        dirV = M(:,3)
        dirH = M(:,2)
      case ("light")
        read(argStr, *, iostat=ios) a, b, c
        light = [a, b, c]
      case ("bounds")
        read(argStr, *, iostat=ios) a, b, c, d, e, f
        bounds = [a, b, c, d, e, f]
      case ("amb")
        read(argStr, *, iostat=ios) a
        a = max(a, ZERO)
        a = min(a, ONE)
        ambient = a

      case ("transp")
        matName = trim(argStr)
        matIdx = mm_matIdx(matName)
        
        if (matIdx == NOT_FOUND) then
          print *,'Material '//trim(matName)//' is not present in the simulation.'
          cycle
        end if

        ! Check mats to see if matIdx is present
        if (any(mats == matIdx)) then
          print *,'Material '//trim(matName)//' is already transparent.'
          cycle
        end if

        n = size(mats)

        allocate(temp(n+1))
        temp(1:n) = mats
        temp(n+1) = matIdx

        call move_alloc(temp, mats)
      
      case ("opaq")
        matName = trim(argStr)
        matIdx = mm_matIdx(matName)
        
        if (matIdx == NOT_FOUND) then
          print *,'Material '//trim(matName)//' is not present in the simulation.'
          cycle
        end if

        ! Check mats to see if matIdx is present
        if (all(mats /= matIdx)) then
          print *,'Material '//trim(matName)//' is already opaque.'
          cycle
        end if

        mats = pack(mats, mats /= matIdx)
      
      case("mats")
        print *,'MATERIAL NAMES:'
        do i = 1, mm_nMat()
          print *,trim(mm_matName(i))
        end do
        cycle

      case("res")
        read(argStr, *, iostat=ios) i, j
    
        if (i < 1 .or. j < 1) then
          print *,'Invalid pixels values provided'
          cycle
        end if

        res(1) = i
        res(2) = j
        deallocate(matIDs, lum)
        allocate(matIDs(res(1), res(2)))
        allocate(lum(res(1), res(2)))

      case ("pos")
        print *,'Camera position: '//numToChar(camera)
        print *,'Centre position: '//numToChar(centre)
        print *,'Light position: '//numToChar(light)
        cycle
      case ("q")
        exit
      case default
        print *,'Unrecognised command'
        cycle
      end select
      
      ! lum contains luminosity values, matIDs identifies which materials were hit
      call self % geom % rayPlot(lum, matIDs, camera, light, M, mats, fov, ambient, bounds)
    
      ! Translate to an image.
      ! Obtain material colours and scale by luminosity
      matIDs = materialColour(matIDs, offset)
      img = matIDs

      ! Scale image by brightness
      img = brightnessScale(img, lum)

      ! Print image
      call imgBmp_toFile(img, outputFile)

      ! Display image
      call execute_command_line("mv "//outputFile//" "//outputLive)

    end do

  end subroutine run

  !!
  !! Initialise from individual components and dictionaries
  !!
  subroutine init(self, dict)
    class(renderPhysicsPackage), intent(inout) :: self
    class(dictionary), intent(inout)        :: dict
    class(dictionary),pointer               :: tempDict, vizDict
    character(nameLen)                      :: geomName
    logical(defBool)                        :: found
    integer(shortInt)                       :: i
    character(nameLen), dimension(:), allocatable :: keys
    character(nameLen)                      :: vizType
    character(100), parameter :: Here ='init (renderPhysicsPackage_class.f90)'

    ! Read the app used to display the live-updated image
    call dict % getOrDefault(self % displayApp, "app", "eog")
    
    ! Build Nuclear Data
    call ndReg_init(dict % getDictPtr("nuclearData"))

    ! Build geometry
    tempDict => dict % getDictPtr('geometry')
    geomName = 'visualGeom'
    call new_geometry(tempDict, geomName)
    self % geomIdx = gr_geomIdx(geomName)
    self % geom    => gr_geomPtr(self % geomIdx)

    ! Find ray dictionary
    if (.not. dict % isPresent('viz')) call fatalError(here,'Must provide viz dict for plotting')
      
    ! Go through dictionary until a subdictionary has a rayPlot
    vizDict => dict % getDictPtr('viz')
    call vizDict % keys(keys)

    found = .false.
    do i = 1, size(keys)
      tempDict => vizDict % getDictPtr(keys(i))

      call tempDict % get(vizType, 'type')
      if (vizType == 'ray') then
        found = .true.
        self % dict = tempDict
      end if
    end do

    if (.not. found) call fatalError(Here,'Must provide a rayPlot dictionary')

  end subroutine init

  !!
  !! Deallocate memory
  !!
  subroutine kill(self)
    class(renderPhysicsPackage), intent(inout) :: self

    ! TODO: This subroutine

  end subroutine kill

end module renderPhysicsPackage_class
