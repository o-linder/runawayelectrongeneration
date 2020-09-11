C######################################################################|
      module calculate_avalanche_growthrate
C----------------------------------------------------------------------|
C     This module is used to calculate the avalanche runaway growth 
C     rate.
C
C     Functions:
C       - Dreicer_growthrate_classic
C       - Dreicer_growthrate_CODE_neural_network
C
C     Function arguments:
C       - Z:    atomic number of each ion species
C       - Z0:   net charge of plasma for each ion species
C       - nj:   density of each ion species (m**-3)
C       - nSpe: number of ion species
C       - T:    electron temperature (eV)
C       - Epar: parallel electric field (V/m)
C
C----------------------------------------------------------------------|
      use double
      use physical_constants, only :
     >  pi
      use Coulomb_logarithms, only :
     >  ln_Lambda_0, ln_Lambda_c
      use collision_frequencies, only :
     >  nu_ee, nuBar_D, nuBar_S
      use electric_fields

      public ::
     >  avalanche_growthrate_classic,
     >  avalanche_growthrate,
     >  p_star

      private ::
     >  dp, E_c, E_ceff_over_Ec, E_D, iter_acc, iter_max, ln_Lambda_0, 
     >  ln_Lambda_c, nu_ee, nuBar_D, nuBar_S, pi

        ! ----- Parameters --------------------------------------------|
      integer, parameter ::
     >  iter_max    = 100               ! maximum number of iterations

      real(kind=dp), parameter ::
     >  iter_acc    = 1.0e-8_dp         ! target accuracy of iterations

      contains

C======================================================================|
      real(kind=dp) function avalanche_growthrate_classic(Z, Z0, nj, 
     >  nSpecies, T, Epar, eps) result (Gamma_av)
C----------------------------------------------------------------------|
C     The flux of secondary electrons to relativistic momenta is 
C     calculated according to Eq. (18) in M.N. Rosenbluth and S.V. 
C     Putvinski, Nucl. Fusion 37, 1355 (1997).
C          
C     Input:
C       Z   : atomic number of each ion species
C       Z0  : net charge of plasma for each ion species
C       nj  : density of each ion species (m**-3)
C       nSpe: number of ion species
C           Z, Z0 and nj must be vectors of length nSpecies
C       T   : temperature (eV)
C       Epar: parallel electric field (V/m)
C       eps : local r/R
C
C     Output:
C       Gamm: flux of secondary runaway electrons (1/s)
C          
C----------------------------------------------------------------------|
        ! ----- Function arguments ------------------------------------|
      integer, intent(in) :: 
     >  nSpecies, Z(nSpecies), Z0(nSpecies)

      real(kind=dp), intent(in) ::
     >  Epar, eps, nj(nSpecies), T

        ! ----- Local variables ---------------------------------------|
      integer ::
     >  i

      real(kind=dp) ::
     >  Ec, gam, ne, nu_ee_th, Zeff

C----------------------------------------------------------------------|
C     Construct required quantities
C----------------------------------------------------------------------|
        ! ----- Free electron density ---------------------------------|
      ne = sum(Z0(:)*nj(:))

        ! ----- Critical electric field -------------------------------|
      Ec = E_c(ne, T)

        ! ----- Thermal electron electron collision frequency ---------|
      nu_ee_th = nu_ee(ne, T)*ln_Lambda_0(ne,T)/ln_Lambda_c(ne, T)
     >  * 7.7431021134304697e-09_dp*T**1.5_dp

        ! ----- Effective charge --------------------------------------|
      Zeff = sum(Z0(:)**2*nj(:))/ne

        ! ----- Geometric factor --------------------------------------|
      gam = 1._dp/(1 + 1.46_dp*sqrt(eps) + 1.72_dp*eps)

C----------------------------------------------------------------------|
C     Calculate growth rate
C----------------------------------------------------------------------|
      Gamma_av = nu_ee_th*sqrt(pi*gam/3._dp/(Zeff+5))*(Epar/Ec-1)
     >  / ln_Lambda_c(ne, T) / sqrt(1 - Ec/Epar + 4*pi*(Zeff+1)**2
     >      / 3._dp/gam/(Zeff+5)/(Epar**2/Ec**2 + 4._dp/gam**2-1))

      return
      end function avalanche_growthrate_classic
C======================================================================|


C======================================================================|
      real(kind=dp) function avalanche_growthrate(Z, Z0, nj, nSpecies, 
     >  T, Epar, B) result (Gamma_av) 
C----------------------------------------------------------------------|
C     The flux of secondary electrons to relativistic momenta is 
C     calculated according to Eq. (14) in L. Hesslow et al., Nucl. 
C     Fusion 59, 084004 (2019).
C          
C     Input:
C       Z   : atomic number of each ion species
C       Z0  : net charge of plasma for each ion species
C       nj  : density of each ion species (m**-3)
C       nSpe: number of ion species
C           Z, Z0 and nj must be vectors of length nSpecies
C       T   : temperature (eV)
C       Epar: parallel electric field (V/m)
C       B   : magnetic field (T)
C
C     Output:
C       Gamm: flux of secondary runaway electrons (1/s)
C          
C----------------------------------------------------------------------|
        ! ----- Function arguments ------------------------------------|
      integer, intent(in) :: 
     >  nSpecies, Z(nSpecies), Z0(nSpecies)

      real(kind=dp), intent(in) ::
     >  B, Epar, nj(nSpecies), T

        ! ----- Local variables ---------------------------------------|
      integer ::
     >  i

      real(kind=dp) ::
     >  ne, ne_tot, pstar

C----------------------------------------------------------------------|
C     Construct required quantities
C----------------------------------------------------------------------|
        ! ----- Free electron density ---------------------------------|
      ne = sum(Z0(:)*nj(:))

        ! ----- Total electron density --------------------------------|
      ne_tot = sum(Z( :)*nj(:)) 

        ! ----- Effective critical momentum for runaway ---------------| 
      pstar = p_star(Z, Z0, nj, nSpecies, T, Epar, B)

C----------------------------------------------------------------------|
C     Calculate growth rate
C----------------------------------------------------------------------|
      Gamma_av  = 5.8667920978923496e+02_dp
     >  * ne_tot/ne * E_c(ne, T) / ln_Lambda_c(ne, T)
     >  * (Epar/E_c(ne, T) - E_ceff_over_E_c(Z, Z0, nj, nSpecies, T,B))
     >  / sqrt(4 + nuBar_D(Z, Z0, nj, nSpecies, T, pstar)
     >           * nuBar_S(Z, Z0, nj, nSpecies, T, pstar)) 

      return
      end function avalanche_growthrate
C======================================================================|


C======================================================================|
      real(kind=dp) function p_star(Z, Z0, nj, nSpecies,
     >  T, Epar, B)
C----------------------------------------------------------------------|
C     Calculates the critical normalised momentum for electron runaway
C     according to L. Hesslow et al., Nucl. Fusion 59, 084004 (2019),
C     just below Eq. (12).
C          
C     Input:
C       Z   : atomic number of each ion species
C       Z0  : net charge of plasma for each ion species
C       nj  : density of each ion species (m**-3)
C       nSpe: number of ion species
C           Z, Z0 and nj must be vectors of length nSpecies
C       T   : temperature (eV)
C       E   : parallel electric field (V/m)
C
C     Output:
C       psta: critical momentum for electron runaway, normalised to me/c
C
C----------------------------------------------------------------------|
        ! ----- Function arguments ------------------------------------|
      integer, intent(in) ::
     >  nSpecies, Z(nSpecies), Z0(nSpecies)

      real(kind=dp), intent(in) ::
     >  B, nj(nSpecies), T, Epar 

        ! ----- Local variables ---------------------------------------|
      integer ::
     >  i

      real(kind=dp) ::
     >  EE_c, ne, p_old

C----------------------------------------------------------------------|
C     Construct required quantities
C----------------------------------------------------------------------|
        ! ----- Free electron density ---------------------------------|
      ne = sum(Z0(:)*nj(:))

        ! ----- Normalised electric field -----------------------------|
      EE_c = max(Epar/E_c(ne, T),
     >           E_ceff_over_E_c(Z, Z0, nj, nSpecies, T, B))

C----------------------------------------------------------------------|
C     Obtain critical momentum iteratively
C----------------------------------------------------------------------|
      p_star = 1._dp

      do i=1,iter_max
        p_old = p_star
        p_star = sqrt(1._dp/EE_c*sqrt(
     >        nuBar_S(Z(:), Z0(:), nj(:), nSpecies, T, p_star)
     >      * nuBar_D(Z(:), Z0(:), nj(:), nSpecies, T, p_star)))

        if (abs(p_star-p_old) .le. iter_acc)
     >      exit 
      enddo

      return
      end function p_star
C======================================================================|

      end module calculate_avalanche_growthrate
C######################################################################|
