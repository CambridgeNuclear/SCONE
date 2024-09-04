
# Add Source Files to the global list
add_sources( ./energyLawENDF_inter.f90
             ./energyLawENDFslot_class.f90
             ./energyLawENDFfactory_func.f90
             ./contTabularEnergy_class.f90
             ./energyEndfPdf/tabularEnergy_class.f90
             ./noEnergy_class.f90
             ./maxwellSpectrum_class.f90
             ./levelScattering_class.f90
             ./evaporationSpectrum_class.f90
             ./multipleEnergyLaws_class.f90
             ./testEnergyLaw_class.f90)

# Add Tests to global list
add_unit_tests( ./Tests/multipleEnergyLaws_test.f90)

