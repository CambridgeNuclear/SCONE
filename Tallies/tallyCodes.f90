!!
!! This module contains codes used in tallies to indentify diffrent events
!!
module tallyCodes

  use numPrecision

  implicit none
  private


  ! List of codes for diffrent reports
  integer(shortInt),parameter,public :: inColl_CODE     = 1000 ,&
                                        outColl_CODE    = 1001 ,&
                                        path_CODE       = 1002 ,&
                                        trans_CODE      = 1003 ,&
                                        hist_CODE       = 1004 ,&
                                        cycleStart_CODE = 1005 ,&
                                        cycleEnd_CODE   = 1006

  ! List of codes for fiffrent particle fates
  integer(shortInt),parameter,public :: abs_FATE  = 5000 ,&
                                        leak_FATE = 5001 ,&
                                        lost_FATE = 5002 ,&
                                        aged_FATE = 5003

end module tallyCodes
