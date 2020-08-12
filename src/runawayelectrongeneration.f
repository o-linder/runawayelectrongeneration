C######################################################################|
      module runawayElectronGeneration
C----------------------------------------------------------------------|
C     This module is used to calculate runaway electron growth rates
C     and related quantities.
C----------------------------------------------------------------------|
      use calculate_Dreicer_growthrate
      use calculate_avalanche_growthrate
      use calculate_electric_fields
      use calculate_collision_frequencies
      use calculate_Coulomb_logarithms

      implicit none

      public  :: 
        ! ----- Runaway electron growth rates -------------------------|
            ! Dreicer growth rates
     >  Dreicer_growthrate_classic,
     >  Dreicer_growthrate_CODE_neural_network,
     >
            ! Avalanche growth rates
     >  avalanche_growthrate_classic,
     >  avalanche_growthrate,
     >
        ! ----- Related quantities ------------------------------------|
            ! Characteristic electric fields
     >  E_ceff_over_E_c,
     >  E_c,
     >  E_D,
     >
            ! Critical momentum
     >  p_star,
     >
            ! Collision frequencies
     >  nuBar_D,
     >  calc_nuBar_D01,
     >  nuBar_S,
     >  calc_nuBar_S01,
     >  nu_ee,
     >
            ! Coulomb logarithms
     >  ln_Lambda_0,
     >  ln_Lambda_c,
     >  ln_Lambda_ee,
     >  ln_Lambda_ei

      contains

      end module runawayElectronGeneration
C######################################################################|
