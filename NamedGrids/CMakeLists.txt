# Add Source Files to the global list
add_sources (./energyGridRegistry_mod.f90
             ./preDefEnergyGrids.f90)

# Add Tests
add_unit_tests(./Tests/energyGridRegistry_test.f90)
