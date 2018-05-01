module dictContent_class

  use numPrecision

  implicit none
  private

  integer(shortInt),parameter,public :: charLen  = max(nameLen,pathLen)
  integer(shortInt),parameter,public :: empty    = 0
  integer(shortInt),parameter,public :: numInt   = 1
  integer(shortInt),parameter,public :: numReal  = 2
  integer(shortInt),parameter,public :: word     = 3
  integer(shortInt),parameter,public :: nestDict = 4
  integer(shortInt),parameter,public :: arrInt   = 5
  integer(shortInt),parameter,public :: arrReal  = 6
  integer(shortInt),parameter,public :: arrWord  = 7



contains


    
end module dictContent_class
