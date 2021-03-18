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
      use collision_frequencies, only :
     >  calc_nuBar_D01, calc_nuBar_S01


      implicit none

      private
      
      public ::
     >  E_ceff_over_E_c,
     >  E_c,
     >  E_D

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
     >  a0, a1, a2, ne, nuBar_D0, nuBar_D1, nuBar_S0, nuBar_S1, 
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
C     Calculate Eceff/Ec analytically
C----------------------------------------------------------------------|
C     With ideas from Vinodh Bandaru (IPP)

        ! ----- Coefficients of cubic equation ------------------------|
        ! a3*x**3 + a2*x**2 + a1*x + a0 = 0
        ! Note, that a3 = 1 in this case and is thus ignored.
      a0 = log(nuBar_D0/(2*nuBar_S1))

      a2 = -2._dp*(nuBar_S0 + nuBar_S1*(1+nuBar_D1/nuBar_D0)*a0)
      a1 = (.5_dp*a2)**2 - (.7_dp+.4_dp*a0)*alpha*nuBar_D0*nuBar_D1 
     >      - nuBar_S1**2
      a0 = -2._dp*tau_rad_inv*nuBar_D0**2

        ! ----- Factors of the solution for x -------------------------|
        ! For the solution, the factors are
        ! a0 -> R = (9*a2*a1 - 27*a0 - 2*a2**3)/54
        ! a1 -> Q = (3*a1 - a2**2)/9
        ! a1 -> D = Q**3 + R**2
      a0 = (9*a2*a1 - 27*a0 - 2*a2**3)/54._dp
      a1 = (3*a1 - a2**2)/9._dp
      a1 = a1**3 + a0**2

        ! ----- Solution for x ----------------------------------------|
        ! The solution is given by:
        ! x = (R + sqrt(D))**(1/3) + (R - sqrt(D))**(1/3) - a2/3
        ! Since it may be that D < 0, complex numbers have to be used
        ! Note, that a1 -> D and a0 -> R
      E_ceff_over_E_c = real(
     >    (a0 + sqrt(cmplx(a1, 0._dp, kind=dp)))**(1._dp/3._dp)
     >  + (a0 - sqrt(cmplx(a1, 0._dp, kind=dp)))**(1._dp/3._dp)
     >  - a2/3._dp)

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
