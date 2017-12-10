module particle_class

  use numPrecision
  use genericProcedures, only : fatalError

  implicit none
  private

  type, public :: particle
  !! Type that describes a particle. Supports both Continious Energy(CE) and Multi-Group (MG) data
  !! but crashes if query is made for CE value of energy for MG particle and vice versa.
    private
    real(kind=defReal),dimension(3) :: r         !! Particle Position
    real(kind=defReal),dimension(3) :: dir       !! Particle Direction
    real(kind=defReal)              :: E         !! Particle Energy
    integer(kind=shortInt)          :: G         !! Particle Energy Group
    real(kind=defReal)              :: w         !! Particle Weight

    integer(kind=shortInt)          :: CurrentCellIndex
    integer(kind=shortInt)          :: CurrentMaterialIndex

    logical(kind=defBool)           :: isDeadFlag
    logical(kind=defBool)           :: isMGFlag
  contains
        !! Public Interface
        ! Constructor
        generic              :: build => buildCE, buildMG
        ! Change Particle State
        generic              :: setEnergy => setEnergyMG, setEnergyCE
        procedure            :: setWeight
        procedure            :: setPosition
        procedure            :: setDirection
        procedure            :: setCell
        procedure            :: setMaterial
        procedure            :: makeMG
        procedure            :: makeCE
        ! Get Data from the particle
        generic              :: Energy => EnergyMG, EnergyCE
        procedure            :: Weight
        procedure            :: Position
        procedure            :: Direction
        procedure            :: Cell
        procedure            :: Material
        procedure            :: isMG
        procedure            :: isDead
        !! Private - Implementation specific procedures
        procedure,private    :: buildCE
        procedure,private    :: buildMG
        procedure,private    :: setEnergyCE
        procedure,private    :: setEnergyMG
        procedure,private    :: EnergyCE
        procedure,private    :: EnergyMG
  end type particle


contains

  function isDead(self) result (isIt)
    class(particle), intent(inout)     :: self
    logical(kind=defBool)              :: isIt
    isIt = self % isDeadFlag
  end function isDead

  function isMG(self) result (isIt)
    class(particle), intent(inout)     :: self
    logical(kind=defBool)              :: isIt
    isIt = self % isMGFlag
  end function isMG

  subroutine Material(self, mat)
    class(particle), intent(inout)       :: self
    integer(kind=shortInt),intent(out)   :: mat
    mat = self % CurrentMaterialIndex
  end subroutine Material

  subroutine Cell(self, cellIndex)
    class(particle), intent(inout)       :: self
    integer(kind=shortInt),intent(out)   :: cellIndex
    cellIndex = self % CurrentCellIndex
  end subroutine Cell

  subroutine Direction(self, dir)
    class(particle), intent(inout)                  :: self
    real(kind=defReal), dimension(3),intent(out)    :: dir
    dir = self % dir
  end subroutine Direction

  subroutine Position(self, r)
    class(particle), intent(inout)                 :: self
    real(kind=defReal), dimension(3),intent(out)   :: r
    r = self % r
  end subroutine Position

  subroutine Weight(self, W)
    class(particle), intent(inout)     :: self
    real(kind=defReal),intent(out)     :: W
    W = self % w
  end subroutine Weight

  subroutine EnergyMG(self, G)
    class(particle), intent(inout)      :: self
    integer(kind=shortInt),intent(out)  :: G

    if (self%isMGFlag) then
      G = self % G
    else
      call fatalError('Energy method of particle (particle_class.f03)', &
                      'Requested group number of CE particle ')
    end if
  end subroutine EnergyMG

  subroutine EnergyCE(self, E)
    class(particle),intent(inout)     :: self
    real(kind=defReal),intent(out)    :: E

    if (self%isMGFlag) then
      call fatalError('Energy method of particle (particle_class.f03)', &
                      'Requested CE energy of MG particle ')
    else
      E = self % E
    end if
  end subroutine EnergyCE


  subroutine setMaterial(self,mat)
    class(particle), intent(inout)     :: self
    integer(kind=shortInt),intent(in)  :: mat
    self % CurrentMaterialIndex = mat
  end subroutine setMaterial

  subroutine setCell(self,cell)
    class(particle), intent(inout)     :: self
    integer(kind=shortInt),intent(in)  :: cell
    self % CurrentCellIndex = cell
  end subroutine setCell

  subroutine setDirection(self,dir)
    class(particle), intent(inout)                   :: self
    real(kind=defReal),dimension(3) ,intent(in)      :: dir
    self % dir = dir
  end subroutine setDirection

  subroutine setPosition(self,r)
    class(particle), intent(inout)                   :: self
    real(kind=defReal),dimension(3) ,intent(in)      :: r
    self % r = r
  end subroutine setPosition

  subroutine setWeight(self,w)
    class(particle), intent(inout)      :: self
    real(kind=defReal), intent(in)      :: w
    self % w = w
  end subroutine setWeight

  subroutine makeMG(self,G)
    class(particle), intent(inout)      :: self
    integer(kind=shortInt), intent(in)  :: G
    self % G = G
    self % isMGFlag = .true.
  end subroutine makeMG

  subroutine makeCE(self,E)
    class(particle), intent(inout)      :: self
    real(kind=defReal),intent(in)       :: E
    self % E = E
    self % isMGFlag = .false.
  end subroutine

  subroutine setEnergyMG(self,G)
    class(particle), intent(inout)     :: self
    integer(kind=shortInt), intent(in) :: G

    if (self%isMGFLag) then
      self % G = G
    else
      call fatalError('setEnergy method of particle (particle_class.f03)', &
                      'An energy of CE particle was set with MG value')
    end if
  end subroutine setEnergyMG

  subroutine setEnergyCE(self,E)
    class(particle), intent(inout)     :: self
    real(kind=defReal), intent(in)     :: E

    if (self%isMGFlag) then
      call fatalError('setEnergy method of particle (particle_class.f03)', &
                      'An energy of MG particle was set with CE value')
    else
      self % E = E
    end if
  end subroutine setEnergyCE

  subroutine buildCE(self,r,dir,E,w,Cell,Mat)
    class(particle), intent(out)                 :: self
    real(kind=defReal),dimension(3),intent(in)   :: r, dir
    real(kind=defReal),intent(in)                :: E, w
    integer(kind=shortInt),intent(in)            :: Cell, Mat

    self % r =r
    self % dir = dir
    self % E = E
    self % w = w
    self % CurrentCellIndex = Cell
    self % CurrentMaterialIndex = Mat

    self % isDeadFlag = .false.
    self % isMGFlag = .false.
  end subroutine
    
  subroutine buildMG(self,r,dir,G,w,Cell,Mat)
    class(particle), intent(out)                 :: self
    real(kind=defReal),dimension(3),intent(in)   :: r, dir
    real(kind=defReal),intent(in)                :: w
    integer(kind=shortInt),intent(in)            :: G
    integer(kind=shortInt),intent(in)            :: Cell, Mat

    self % r =r
    self % dir = dir
    self % G = G
    self % w = w
    self % CurrentCellIndex = Cell
    self % CurrentMaterialIndex = Mat

    self % isDeadFlag = .false.
    self % isMGFlag = .true.
  end subroutine

end module particle_class
