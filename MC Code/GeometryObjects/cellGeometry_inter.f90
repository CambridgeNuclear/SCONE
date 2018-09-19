module cellGeometry_inter

  use numPrecision
  use coord_class,    only : coordList
  use geometry_inter, only : geometry

  implicit none
  private

  !!
  !! Simplified Cell Geometry Model
  !!   -> Geometry is composed of cells, boundary and a single trip surface
  !!   -> Cells can contain a single material or a number of other cells
  !!   -> Each cell is assigned with an unique ID
  !!   -> Geometry spans the entire space, which is divided into
  !!      domain and outside indicated by a special material
  !!   -> Boundary is the interface surface between domain and outside
  !!   -> Boundary is a single surface
  !!   -> Boundary conditions are handeled by coordinate transformations
  !!   -> Vacuum, reflective and periodic BC are supported
  !!   -> Calculation domain must be convex
  !!   -> Geometry can contain a single trip surface
  !!   -> Intersection of trip surface with boundary cannot be a surface
  !!
  type, public,extends(geometry),abstract :: cellGeometry
    private
  contains
    procedure(move),deferred       :: move
    procedure(teleport),deferred   :: teleport
    procedure(moveGlobal),deferred :: moveGlobal
  end type cellGeometry

  abstract interface

    !!
    !! Given coordinates placed in a geometry move point through the geometry
    !! Move by maxDist. Stop at the boundaries between unique cell IDs or at the trip surface
    !!
    !! If during motion boundary is crossed transformations are applied and particle is stopped
    !! If particle is stopped at the trip surface tripFlag is set to .true.
    !! If particle stoped at the collision isColl is set to .true.
    !!
    !! maxDist is set to (maxDist - "distance moved") at the end of execution
    !!
    !! NOTE:
    !!  -> if coords is not placed in the geometry behaviour is unspecified
    !!  -> if maxDist < 0.0 behaviour is unspecified
    !!
    subroutine move(self,coords,maxDist,isColl,tripFlag)
      import :: cellGeometry, &
                coordList,&
                defReal, &
                defBool
      class(cellGeometry), intent(inout)    :: self
      type(coordList), intent(inout)        :: coords
      real(defReal),intent(inout)           :: maxDist
      logical(defBool), intent(out)         :: isColl
      logical(defBool),optional,intent(out) :: tripFlag
    end subroutine move

    !!
    !! Given coordinates transport point outside the geometry by distance dist
    !!
    !! If during motion boundary is crossed transformation are applied
    !! Particle is NOT stoped but transport continues until total
    !! transport length is equal to dist.
    !!
    !! NOTE:
    !!  -> if coords is not placed in the geometry behaviour is unspecified
    !!  -> if maxDist < 0.0 behaviour is unspecified
    !!
    subroutine teleport(self,coords,dist)
      import :: cellGeometry, &
                coordList, &
                defReal
      class(cellGeometry), intent(inout) :: self
      type(coordList), intent(inout)     :: coords
      real(defReal), intent(in)          :: dist
    end subroutine teleport


    !!
    !! Given coordinates move point in global geometry level
    !! Move by maxDist. Stop at the boundary or at the trip surface
    !!
    !! If during motion boundary is crossed transformations are applied and particle is stopped
    !! If particle is stopped at the trip surface tripFlag is set to .true.
    !! maxDist is set to (maxDist - "distance moved") at the end of execution
    !!
    !! NOTE:
    !!  -> if maxDist < 0.0 behaviour is unspecified
    !!
    subroutine moveGlobal(self,coords,maxDist,tripFlag)
      import :: cellGeometry, &
                coordList, &
                defReal, &
                defBool
      class(cellGeometry), intent(inout)    :: self
      type(coordList), intent(inout)        :: coords
      real(defReal), intent(inout)          :: maxDist
      logical(defBool),optional,intent(out) :: tripFlag
    end subroutine moveGlobal

  end interface

    
end module cellGeometry_inter
