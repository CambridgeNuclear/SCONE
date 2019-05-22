module reactionHandle_inter

  implicit none
  private

  !!
  !! Abstract reaction handle
  !!
  !! Does nothing or has nothing by itself. It only exists to group all types or
  !! reaction interfaces into a single object fimaly
  !!
  type, public,abstract :: reactionHandle
    private
  end type reactionHandle

contains
    
end module reactionHandle_inter
