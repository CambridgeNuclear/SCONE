program test

  use numPrecision
  use RNG_class,               only : RNG
  use aceCard_class,           only : aceCard
  use ceNeutronDatabase_inter, only : ceNeutronDatabase
  use aceNeutronNuclide_class, only : aceNeutronNuclide
  use neutronXsPackages_class, only : neutronMicroXSs

  implicit none
  type(aceNeutronNuclide), target   :: nuc
  type(aceCard)                     :: ACE
  class(ceNeutronDatabase), pointer :: database => null()
  type(neutronMicroXSs) :: xss
  type(RNG) :: R
  integer(shortInt) :: idx
  real(defReal) :: f

  ! Build ACE library
  call ACE % readFromFile('./IntegrationTestFiles/92233JEF311.ace', 1)

  ! Build nuclide
  call nuc % init(ACE, 1, database)

  print *, size(nuc % MTData)
  call nuc % search(idx, f, 5.0_defReal)
  call nuc % microXSs(xss, idx, f)

  print *, xss % total, xss % inelasticScatter, xss % fission

end program test




