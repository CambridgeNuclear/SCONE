module couplingAdmin_iTest

  use numPrecision
  use universalVariables
  use dictionary_class,    only : dictionary
  use dictParser_func,     only : fileToDict
  use fieldFactory_func,   only : new_field
  use couplingAdmin_class, only : couplingAdmin
  use funit

  implicit none

contains

  !!
  !! Coupling admin integration test
  !!
@Test
  subroutine test_coupling()
    type(couplingAdmin)           :: coupler
    character(pathLen), parameter :: path = './IntegrationTestFiles/Coupling/test_coupling'
    character(pathLen), parameter :: sendPath = './IntegrationTestFiles/Coupling/toDriver.dat'
    character(pathLen), parameter :: recvPath = './IntegrationTestFiles/Coupling/toSCONE.dat'
    character(pathLen), parameter :: tallyPath = './IntegrationTestFiles/Coupling/sconeOut.m'
    type(dictionary)              :: dict
    integer(shortInt)             :: it, unit
    logical(defBool)              :: exists, correctFile
    character(nameLen)            :: line

    ! Create a recv file to prevent SCONE from hanging
    open(newunit=unit, file=recvPath, status="replace", action="write")
    write(unit, '(A)') 'SIGUSR1'
    close(unit)

    ! Attempt a coupling update before initialisation: no files should be produced
    ! or deleted
    it = 5
    call coupler % couple(it)

    inquire(file=recvPath, exist=exists)
    @assertTrue(exists)
    inquire(file=tallyPath, exist=exists)
    @assertFalse(exists)
    inquire(file=sendPath, exist=exists)
    @assertFalse(exists)

    ! Initialise
    call fileToDict(dict, path)
    call coupler % init(dict)

    ! Feed in different iteration numbers and see whether files are produced
    ! This one should not result in any information exchange
    it = 1
    call coupler % couple(it)

    inquire(file=recvPath, exist=exists)
    @assertTrue(exists)
    inquire(file=tallyPath, exist=exists)
    @assertFalse(exists)
    inquire(file=sendPath, exist=exists)
    @assertFalse(exists)

    ! This one should be coupled
    it = 5
    call coupler % couple(it)

    inquire(file=recvPath, exist=exists)
    @assertFalse(exists)
    inquire(file=tallyPath, exist=exists)
    @assertTrue(exists)
    inquire(file=sendPath, exist=exists)
    @assertTrue(exists)

    ! Send file should contain SIGUSR1
    open(newunit=unit, file=sendPath, status='old', action='read')
    read(unit,'(A)') line
    close(unit)
    correctFile = (trim(line) == 'SIGUSR1')
    @assertTrue(correctFile)

    ! Reset the files
    open(newunit=unit, file=tallyPath, status='old')
    close(unit, status='delete')
    open(newunit=unit, file=sendPath, status='old')
    close(unit, status='delete')

    ! This one should end coupling due to hitting the maximum iteration
    ! No recv file is necessary
    it = 10
    call coupler % couple(it)

    inquire(file=recvPath, exist=exists)
    @assertFalse(exists)
    inquire(file=tallyPath, exist=exists)
    @assertTrue(exists)
    inquire(file=sendPath, exist=exists)
    @assertTrue(exists)

    ! Send file should contain SIGTERM
    open(newunit=unit, file=sendPath, status='old', action='read')
    read(unit,'(A)') line
    close(unit)
    correctFile = (trim(line) == 'SIGTERM')
    @assertTrue(correctFile)

    ! Reset the files
    open(newunit=unit, file=tallyPath, status='old')
    close(unit, status='delete')
    open(newunit=unit, file=sendPath, status='old')
    close(unit, status='delete')
    open(newunit=unit, file=recvPath, status="replace", action="write")
    write(unit, '(A)') 'SIGTERM'
    close(unit)

    ! Coupling should not occur
    it = 15
    call coupler % couple(it)

    inquire(file=recvPath, exist=exists)
    @assertTrue(exists)
    inquire(file=tallyPath, exist=exists)
    @assertFalse(exists)
    inquire(file=sendPath, exist=exists)
    @assertFalse(exists)

    ! Reinitialise the coupler, ensuring that no coupling occurs successfully
    ! until reinitialisation
    call coupler % kill()
    it = 5
    call coupler % couple(it)

    inquire(file=recvPath, exist=exists)
    @assertTrue(exists)
    inquire(file=tallyPath, exist=exists)
    @assertFalse(exists)
    inquire(file=sendPath, exist=exists)
    @assertFalse(exists)
    
    call coupler % init(dict)

    ! Given the recv file, coupling should occur on this iteration, but not on the 10th
    it = 5
    call coupler % couple(it)
    
    inquire(file=recvPath, exist=exists)
    @assertFalse(exists)
    inquire(file=tallyPath, exist=exists)
    @assertTrue(exists)
    inquire(file=sendPath, exist=exists)
    @assertTrue(exists)
    
    ! Send file should contain SIGUSR1
    open(newunit=unit, file=sendPath, status='old', action='read')
    read(unit,'(A)') line
    close(unit)
    correctFile = (trim(line) == 'SIGUSR1')
    @assertTrue(correctFile)
    
    ! Reset the files - no need for the recv file
    open(newunit=unit, file=tallyPath, status='old')
    close(unit, status='delete')
    open(newunit=unit, file=sendPath, status='old')
    close(unit, status='delete')

    it = 10
    call coupler % couple(it)
    
    inquire(file=recvPath, exist=exists)
    @assertFalse(exists)
    inquire(file=tallyPath, exist=exists)
    @assertFalse(exists)
    inquire(file=sendPath, exist=exists)
    @assertFalse(exists)

    ! Kill coupling admin
    call coupler % kill()

  end subroutine test_coupling

end module couplingAdmin_iTest
