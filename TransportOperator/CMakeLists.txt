# Add Source Files to the global list
add_sources(./transportOperator_inter.f90
            ./transportOperatorFactory_func.f90
            ./transportOperatorDT_class.f90
          #  ./transportOperatorDynamicDT_class.f90
            ./transportOperatorST_class.f90
            ./transportOperatorHT_class.f90
            ./transportOperatorTime_class.f90
            ./transportOperatorTimeHT_class.f90
            ./Grid/trackingGrid_class.f90)
