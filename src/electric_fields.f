C######################################################################|
      module electric_fields
C----------------------------------------------------------------------|
C     This module is used to calculate characteristic electric fields.
C   
C     Functions:
C       - E_ceff_ocer_E_c(.):       normalised effective critical field
C       - E_c(ne, T):               critical electric field
C       - E_D(ne, T):               Dreicer electric field
C
C     Function arguments:
C       - Z:    atomic number of each ion species
C       - Z0:   net charge of plasma for each ion species
C       - nj:   density of each ion species (m**-3)
C       - nSpe: number of ion species
C       - T:    electron temperature (eV)
C       - B:    magnetic field (T)
C       - ne:   free electron density (m**-3)
C
C----------------------------------------------------------------------|
      use double
      use physical_constants, only :
     >  alpha
      use Coulomb_logarithms, only :
     >  ln_Lambda_0, ln_Lambda_c
      use collision_frequencies

      implicit none
      
      public ::
     >  E_ceff_over_E_c,
     >  E_c,
     >  E_D

      private ::
     >  alpha, dp, iter_acc, iter_max, ln_Lambda_0, ln_Lambda_c

        ! ----- Parameters --------------------------------------------|
      integer, parameter ::
     >  iter_max    = 100               ! maximum number of iterations

      real(kind=dp), parameter ::
     >  iter_acc    = 1.0e-8_dp         ! target accuracy of iterations

      contains

C======================================================================|
      real(kind=dp) function E_ceff_over_E_c(Z, Z0, nj, nSpecies,
     >  T, B)
C----------------------------------------------------------------------|
C     The effective critical field is calculated as in 
C     L. Hesslow et al., Plasma Phys. Control. Fusion 60, 074010 (2018) 
C     in Eq. (23)-(24) for a plasma composition with any fully ionised 
C     species in combination with partially ionized neon and/or argon.
C
C       Input:
C           Z   : atomic number of each ion species
C           Z0  : net charge of plasma for each ion species
C           nj  : density of each ion species (m**-3)
C           nSpe: number of ion species
C               Z, Z0 and nj must be vectors of length nSpecies
C           T   : temperature (eV)
C           B   : magnetic field (T)
C
C       Output:
C           Eceff/Ec, where Ec = ne e^3 lnLc/(4 pi epsilon_0^2 m_e c^2)
C
C----------------------------------------------------------------------|
        ! ----- Function arguments ------------------------------------|
      integer, intent(in) ::
     >  nSpecies, Z(nSpecies), Z0(nSpecies)

      real(kind=dp), intent(in) ::
     >  B, nj(nSpecies), T

        ! ----- Local variables ---------------------------------------|
      integer ::
     >  i

      real(kind=dp) ::
     >  e0, e1, e2, E_old, ne, nuBar_D0, nuBar_D1, nuBar_S0, nuBar_S1, 
     >  tau_rad_inv

C----------------------------------------------------------------------|
C     Construct required quantities
C----------------------------------------------------------------------|
        ! ----- Free electron density ---------------------------------|
      ne = sum(Z0(:)*nj(:))
      
        ! ----- Generalised deflection frequency ----------------------|
      call calc_nuBar_D01(Z, Z0, nj, nSpecies, T, nuBar_D0, nuBar_D1)

        ! ----- Slowing down frequency --------------------------------|
      call calc_nuBar_S01(Z, Z0, nj, nSpecies, T, nuBar_S0, nuBar_S1)
  
        ! ----- Synchrotron radiation term ----------------------------|
      tau_rad_inv = B**2/(15.44_dp*ln_Lambda_c(ne, T))/(1.e-20_dp*ne)

C----------------------------------------------------------------------|
C     Calculate normalised Eceff iteratively
C----------------------------------------------------------------------|
        ! ----- Coefficients ------------------------------------------|
      e2 = log(nuBar_D0/(2*nuBar_S1))

      e0 = nuBar_S0 + nuBar_S1*(1+nuBar_D1/nuBar_D0)*e2
      e1 = (.7_dp+.4_dp*e2)*alpha*nuBar_D0*nuBar_D1 + nuBar_S1**2
      e2 = 2*tau_rad_inv*nuBar_D0**2

        ! ----- Initial guess -----------------------------------------|
      E_ceff_over_E_c  = ln_Lambda_c(ne, T)/ln_Lambda_0(ne, T)
     >  * sum(Z(:)*nj(:))/ne

        ! ----- Iterative calculation ---------------------------------|
      do i=1,iter_max
        E_old = E_ceff_over_E_c
        E_ceff_over_E_c = e0 + sqrt(e1 + e2/E_old)

        if (abs(E_ceff_over_E_c-E_old) .le. iter_acc) 
     >      exit 
      enddo

      return
      end function E_ceff_over_E_c
C======================================================================|


C======================================================================|
      real(kind=dp) function E_c(ne, T)
C----------------------------------------------------------------------|
C     Calculates the critical electric field E_c:
C       E_c = me*c/e/tauc = ne*e**3*ln Lambda_c/(4*pi*eps0**2*me*c**2).
C
C     Input:
C       ne  : free electron density (m**-3)
C       T   : electron temperature (eV)
C
C     Output:
C       E_c : critical electric field (V)
C
C----------------------------------------------------------------------|
        ! ----- Function arguments ------------------------------------|
      real(kind=dp), intent(in) ::
     >  ne, T

C----------------------------------------------------------------------|
      E_c = 5.09909908817339515e-23_dp*ne*ln_Lambda_c(ne, T)

      return
      end function E_c
C======================================================================|


C======================================================================|
      real(kind=dp) function E_D(ne, T)
C----------------------------------------------------------------------|
C     Calculates the Dreicer electric field E_D:
C       E_D = ne*e**3*ln Lambda/(4*pi*eps0**2*T).
C
C     Input:
C       ne  : free electron density (m**-3)
C       T   : electron temperature (eV)
C
C     Output:
C       E_D : Dreicer electric field (V)
C
C----------------------------------------------------------------------|
        ! ----- Function arguments ------------------------------------|
      real(kind=dp), intent(in) ::
     >  ne, T

C----------------------------------------------------------------------|
      E_D = 2.6056342609758366e-17_dp*ne*ln_Lambda_0(ne, T)/T

      return
      end function E_D
C======================================================================|


      end module electric_fields
C######################################################################|
