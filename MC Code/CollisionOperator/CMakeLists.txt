# Add Source Files to the global list 
add_sources( #./collisionOperator_class.f90 
             ./scatteringKernels_func.f90
             #./MGOperators/collisionOperatorMG_class.f90
             ./collisionOperatorBase_inter.f90
             ./perNuclideCollisionOpCE_class.f90
             ./perNuclideImplicitCaptureCE_class.f90
             ./perMaterialCollisionOpMG_class.f90)
