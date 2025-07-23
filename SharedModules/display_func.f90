module display_func
  use numPrecision
  use universalVariables, only : MAX_COL
  use iso_fortran_env,    only : compiler_version
  use mpi_func,           only : isMPIMaster, getMPIWorldSize
  use openmp_func,        only : ompGetMaxThreads
  implicit none

contains

  !!
  !! Prints Scone ACII Header
  !!
  subroutine printStart()
    if (isMPIMaster()) then
      print *, repeat(" ><((((*> ", MAX_COL / 10)
      print *, ''
      print * ,"        _____ __________  _   ________  "
      print * ,"       / ___// ____/ __ \/ | / / ____/  "
      print * ,"       \__ \/ /   / / / /  |/ / __/     "
      print * ,"      ___/ / /___/ /_/ / /|  / /___     "
      print * ,"     /____/\____/\____/_/ |_/_____/     "
      print * , ''
      print * , ''
      print * , "Compiler Info :   ", compiler_version()
#ifdef _OPENMP
      print '(A, I4)', " OpenMP Threads: ", ompGetMaxThreads()
#endif
#ifdef MPI
      print '(A, I4)', " MPI Processes:  ", getMPIWorldSize()
#endif
      print *
      print *, repeat(" <*((((>< ", MAX_COL / 10)
    endif
    ! TODO: Add extra info like date & time

  end subroutine printStart


  !!
  !! Prints line of fishes swimming right with an offset
  !!
  !! Args:
  !!   offset [in] : Number of characters to shift the line right
  !!
  subroutine printFishLineR(offset)
    integer(shortInt),intent(in) :: offset
    integer(shortInt)            :: offset_L
    character(MAX_COL), dimension(10), parameter :: lines = [&
    "" // repeat(" ><((((*> ", MAX_COL / 10 - 1) // " ><((((*> ", &
    " " // repeat(" ><((((*> ", MAX_COL / 10 - 1) // " ><((((*>",&
    "> " // repeat(" ><((((*> ", MAX_COL / 10 - 1) // " ><((((*",&
    "*> " // repeat(" ><((((*> ", MAX_COL / 10 - 1) // " ><((((",&
    "(*> " // repeat(" ><((((*> ", MAX_COL / 10 - 1) // " ><(((",&
    "((*> " // repeat(" ><((((*> ", MAX_COL / 10 - 1) // " ><((",&
    "(((*> " // repeat(" ><((((*> ", MAX_COL / 10 - 1) // " ><(",&
    "((((*> " // repeat(" ><((((*> ", MAX_COL / 10 - 1) // " ><",&
    "<((((*> " // repeat(" ><((((*> ", MAX_COL / 10 - 1) // " >",&
    "><((((*> " // repeat(" ><((((*> ", MAX_COL / 10 - 1) // " "]

    offset_L = modulo(offset, 10)

    if (isMPIMaster()) then
      print *, lines(offset_L+1)
    endif

  end subroutine  printFishLineR

  !!
  !! Prints a message to the screen
  !!
  !! Needs to be used in place of `print *` to ensure that the message is only
  !! printed by the master process when using MPI
  !!
  subroutine statusMsg(msg)
    character(*), intent(in) :: msg
    if (isMPIMaster()) then
      print *, msg
    endif
  end subroutine statusMsg

  !!
  !! Print a section start header
  !!
  !!  " /\/\ Section Name /\/\"
  !!
  !! Args:
  !!   name [in] : Name of the section
  !!
  subroutine printSectionStart(name)
    character(*), intent(in) :: name

    if (isMPIMaster()) then
      print *, "/\/\ " // name // " /\/\"
    endif

  end subroutine printSectionStart

  !!
  !! Print a section end header
  !!
  !!  " \/\/ Section Name \/\/"
  !!
  !! Args:
  !!   name [in] : Name of the section
  !!
  subroutine printSectionEnd(name)
    character(*), intent(in) :: name

    if (isMPIMaster()) then
      print *, "\/\/ " // name // " \/\/"
    endif

  end subroutine printSectionEnd

  !!
  !! Prints a separator line
  !!
  !! "<><><><><>..." to max column width
  !!
  subroutine printSeparatorLine()
    if (isMPIMaster()) then
      print *, repeat("<>", MAX_COL / 2)
    endif
  end subroutine printSeparatorLine

end module display_func
