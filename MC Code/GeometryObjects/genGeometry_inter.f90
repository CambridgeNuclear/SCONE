module genGeometry_inter

  use numPrecision
  use geometry_inter, only : geometry

  implicit none
  private

  !! *** Specification requires further revision and clear interface definition ***
  !!
  !! Flexible Gen(eral) Geometry Model
  !!   -> Geometry is composed of cells, surfaces and universes
  !!   -> Universe divides entire space into cells
  !!   -> Cells can contain a single material or an universe
  !!   -> Every cell is assigned with an uniqe ID
  !!   -> Base universe divides space into domain and outside indicated by a special material
  !!   -> Two cells are neighbours when ther intersection is a surface.
  !!   -> Each cell is a union of closed surface halfspaces. It is allowed that neighbouring cells
  !!      do not share a surface in their definition!
  !!   -> Particle at the surface is in a half-space it moves into (based on its direction)
  !!   -> Particle at a plane with direction parallel to the plane is lost
  !!   -> Surfaces can be simple or compound.
  !!   -> Compound surfaces are a collection of simple surfaces
  !!   -> Trip surfaces can be present in base universe only. They cannot be used to define cells.
  !!   -> Simple surface, that is not part of a compound surface, can have a single surface operator
  !!   -> Compound surface can have a single surface operator
  !!   -> Operators act on a particle crossing a surface
  !!   -> Multiple surfaces can be crossed at a single point, but only one operator is applied
  !!   -> Surface crossings are tallied before application of the operator
  !!   -> Operators can change everything in a public interface of a particle
  !!   -> At the end of an operator action, particle must be placed in its current geometry
  !!   -> On multiple surface crossings operator of highest priority surface is applied. If this
  !!      surface has no operator no operators are applied.
  !!   -> Highest priority surface is a trip surface or the highest level surface used to
  !!      define current particle cell if no trip surface is crossed.
  !!   -> On crossing two highest priority surfaces (corner, overlapping section), it is undefined
  !!      which operator will be applied
  !!   -> If after application of an operator particle is outside geometry it is lost.
  !!

  type,public,extends(geometry),abstract :: genGeometry

  contains
    !procedure(move), deferred       :: move
    !procedure(moveGlobal), deferred :: moveGlobal
  end type genGeometry

contains


end module genGeometry_inter
