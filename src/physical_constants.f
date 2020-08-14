C######################################################################|
      module physical_constants
C----------------------------------------------------------------------|
C     Provides several, mostyl physical constants.
C----------------------------------------------------------------------|
      use double

      implicit none

      public ::
     >  e,
     >  epsilon_0,
     >  alpha, 
     >  
     >  m_e, 
     >  m_ec, 
     >
     >  pi

      private ::
     >  dp

        ! ----- Parameters --------------------------------------------|
      real(kind=dp), parameter ::
     >  e       = 1.602176634e-19_dp,       ! electron charge (C)
     >  epsilon_0 = 8.8541878128e-12_dp,    ! vacuum permittivity (F/m)
     >  alpha   = 7.2973525693e-3_dp,       ! fine-structure constant
     >
     >  m_e     = 9.1093837015e-31_dp,      ! electron mass (kg)
     >  m_ec    = 5.1099895000e+5_dp,       ! electron mass (eV)
     >
     >  pi      = 3.1415926535897931_dp     ! this one you should know

      contains

      end module physical_constants
C######################################################################|
