program test
  use numPrecision

  use nuclearData_inter,       only : nuclearData
  use IOdictionary_class,      only : IOdictionary
  use dictionary_class,        only : dictionary
  use nuclearDataRegistry_mod, only : build_nuclearData, getHandlePtr
  use cellGeometry_inter,      only : cellGeometry
  use geometryFactory_func,    only : new_cellGeometry_ptr
  implicit none

  type(IOdictionary)                          :: inDict
  class(dictionary),pointer                   :: dict
  character(nameLen),dimension(:),allocatable :: matNames
  character(nameLen)                          :: name1, name2
  class(nuclearData), pointer                 :: nucDat
  class(cellGeometry), pointer                :: geom

  
  ! Load input dictionary from a file 
  call inDict % initFrom('./InputFiles/Jakob_input')
  
  ! Retrieve sub-dictionary defining materials & nuclear data 
  dict => inDict % getDictPtr('nuclearData')
  
  ! Build all materials for all nuclear data representations  
  call build_nuclearData(dict)

  name1 = 'uo2'
  name2 = 'water'

  name1 = 'ce2'
  ! Get pointer to defined nuclearData "ce2" of type "datalessMaterials" 
  nucDat => getHandlePtr(name1)
  
  ! Get pointer to sub-dictionary defining geometry and build geometry 
  dict => inDict % getDictPtr('geometry')
  geom => new_cellGeometry_ptr(dict, nucDat)

end program test


