module geometry_inter

  use numPrecision
  use coord_class, only : coordList

  implicit none
  private

  !!
  !! Base class for all geometry interfaces
  !! Specifies a basic user interface for common geometry operations
  !!
  type, public,abstract :: geometry
    private
  contains
    procedure(placeCoord), deferred :: placeCoord
    procedure(whatIsAt), deferred   :: whatIsAt
    procedure(bounds), deferred     :: bounds
    procedure(slicePlot), deferred  :: slicePlot
    procedure(voxelPlot), deferred  :: voxelPlot
  end type geometry

  abstract interface

    !!
    !! Places coordinate list into geometry
    !! Finds unique cell and material as well as coordinates at all intermediate levels
    !! Returns an error if coordList is uninitialised
    !!
    subroutine placeCoord(self,coords)
      import :: geometry,&
                coordList
      class(geometry), intent(in)    :: self
      type(coordList), intent(inout) :: coords
    end subroutine placeCoord

    !!
    !! Given position in a geometry
    !! Return material index and unique cell ID under point r
    !! NOTE: Implementation can choose how to deal with points exactly at the surfaces
    !!
    subroutine whatIsAt(self,r,matIdx,uniqueID)
      import :: geometry, &
                defReal,&
                shortInt
      class(geometry), intent(in)            :: self
      real(defReal),dimension(3),intent(in)  :: r
      integer(shortInt),intent(out)          :: matIdx
      integer(shortInt), intent(out)         :: uniqueID
    end subroutine whatIsAt

    !!
    !! Return bounds of the geometry domain in cartesian co-ordinates
    !! bounds = [x_min, x_max, y_min, y_max, z_min, z_max]
    !! If geometry is infinate in a given axis direction * then:
    !! *_min = *_max = ZERO
    !!
    function bounds(self)
      import :: geometry,&
                defReal
      class(geometry), intent(in) :: self
      real(defReal),dimension(6)  :: bounds
    end function bounds

    !!
    !! Produce a 2D plot of a geometry
    !! Resolution is determined by a size of input matrix colorMat
    !! By default extent of a plot is determined by bounds of the domain and offset is [0,0,0]
    !!
    !! NOTES:
    !! -> what = {"material","cellID"} determines if matIdx is put in colorMat or unique cellID
    !! -> dir  = {"x","y","z"} specifies direction of the plot
    !! -> centre allows to set offset of a plane
    !! -> width sets well... width in each direction of the plane width(1) is either x or y
    !!
    subroutine slicePlot(self, colorMat, centre, dir, what, width)
      import :: geometry, &
                shortInt,&
                defReal
      class(geometry),intent(in)                       :: self
      integer(shortInt),dimension(:,:), intent(inout)  :: colorMat
      real(defReal), dimension(3), intent(in)          :: centre
      character(1), intent(in)                         :: dir
      character(*), intent(in)                         :: what
      real(defReal), dimension(2), optional,intent(in) :: width
    end subroutine slicePlot

    !!
    !! Produce a voxel 3D plot of geometry
    !! Resolution is determined by a size of input voxelMat
    !! By default bounds of the voxel plot correspond to bounds of geometry
    !!
    !! NOTES:
    !! -> what = {"material","cellID"} determines if matIdx is put in voxel Mat or unique cellID
    !! -> centere and optional width specify extent of plot
    !! -> Voxel plot is always a box and it is axis aligned
    !!
    subroutine voxelPlot(self,voxelMat,what,center,width)
      import :: geometry, &
                shortInt, &
                defReal
      class(geometry),intent(in)                       :: self
      integer(shortInt),dimension(:,:,:),intent(inout) :: voxelMat
      character(*),intent(in)                          :: what
      real(defReal),dimension(3),intent(in)            :: center
      real(defReal),dimension(3),optional,intent(in)   :: width
   end subroutine voxelPlot

  end interface

end module geometry_inter
