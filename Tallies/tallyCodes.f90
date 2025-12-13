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
                                        spawn_CODE      = 1004 ,&
                                        hist_CODE       = 1005 ,&
                                        cycleStart_CODE = 1006 ,&
                                        cycleEnd_CODE   = 1007

  ! List of codes for diffrent particle fates
  integer(shortInt),parameter,public :: no_FATE   = 5000 ,&
                                        abs_FATE  = 5001 ,&
                                        leak_FATE = 5002 ,&
                                        lost_FATE = 5003 ,&
                                        aged_FATE = 5004

end module tallyCodes
