C######################################################################|
      module calculate_Dreicer_growthrate
C----------------------------------------------------------------------|
C     This module is used to calculate the Dreicer runaway growth rate.
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
      use collision_frequencies, only :
     >  nu_ee
      use electric_fields, only :
     >  E_c, E_D

      implicit none

      public ::
     >  Dreicer_growthrate_classic,
     >  Dreicer_growthrate_CODE_neural_network

      private ::
     >  dp, E_c, E_D, nu_ee, pi

        ! ----- Parameters --------------------------------------------|
      real(kind=dp), parameter ::
     >  pi = 3.1415926535897931_dp

      contains

C======================================================================|
      real(kind=dp) function Dreicer_growthrate_classic(Z, Z0, nj, 
     >  nSpecies, T, Epar) result (Gamma_D) 
C----------------------------------------------------------------------|
C     The flux of primary electrons to relativistic momenta is 
C     calculated according to Eqs. (62)-(64) in
C     J.W. Connor and R.J. Hastie, Nucl. Fusion 15, 415 (1975).
C          
C     Input:
C       Z   : atomic number of each ion species
C       Z0  : net charge of plasma for each ion species
C       nj  : density of each ion species (m**-3)
C       nSpe: number of ion species
C           Z, Z0 and nj must be vectors of length nSpecies
C       T   : temperature (eV)
C       Epar: parallel electric field (V/m)
C
C     Output:
C       Gamm: flux of primary runaway electrons (1/s)
C          
C----------------------------------------------------------------------|
        ! ----- Parameters --------------------------------------------|
      real(kind=dp), parameter ::
     >  CR = 1.0_dp

        ! ----- Function arguments ------------------------------------|
      integer, intent(in) :: 
     >  nSpecies, Z(nSpecies), Z0(nSpecies)

      real(kind=dp), intent(in) ::
     >  Epar, nj(nSpecies), T

        ! ----- Local variables ---------------------------------------|
      integer ::
     >  i

      real(kind=dp) ::
     >  EE_c, EE_D, gam, h, lambda, ne, Zeff

C----------------------------------------------------------------------|
C     Construct required quantities
C----------------------------------------------------------------------|
        ! ----- Free electron density ---------------------------------|
      ne = sum(Z0(:)*nj(:))

        ! ----- Critical electric field -------------------------------|
      EE_c = Epar/E_c(ne, T)

        ! ----- Dreicer electric field --------------------------------|
      EE_D = Epar/E_D(ne, T)

        ! ----- Effective charge --------------------------------------|
      Zeff = sum(Z0(:)**2*nj(:))/ne

C----------------------------------------------------------------------|
C     Calculate growth rate
C----------------------------------------------------------------------|
      Gamma_D = 0._dp
      if (EE_c .le. 1._dp) return

      h = (EE_c - (Zeff-7)/(Zeff+1) + 2*sqrt(EE_c/(EE_c-1))
     >  *(EE_c-2))*(Zeff+1)/(16*(EE_c-1))

      lambda = 8*EE_c**2*(1 - 1._dp/(2*EE_c) - sqrt(1 - 1._dp/EE_c))

      gam = EE_c/2*sqrt((Zeff+1)/(EE_c-1))*(pi/2 - asin(1-2._dp/EE_c))

      Gamma_D = CR*nu_ee(ne, T)*EE_D**(-h)
     >  *exp(-lambda/(4*EE_D) - gam/sqrt(EE_D))

      return
      end function Dreicer_growthrate_classic
C======================================================================|


C======================================================================|
      real(kind=dp) function Dreicer_growthrate_CODE_neural_network(
     >   Z, Z0, nj, nSpecies, T, Epar) result (Gamma_D)
C----------------------------------------------------------------------|
C     The flux of primary electrons to relativistic momenta is 
C     calculated according to Eq. (3.3) in L. Hesslow et al., J. Plasma
C     Phys. 85, 475850601 (2019).
C          
C     Input:
C       Z   : atomic number of each ion species
C       Z0  : net charge of plasma for each ion species
C       nj  : density of each ion species (m**-3)
C       nSpe: number of ion species
C           Z, Z0 and nj must be vectors of length nSpecies
C       T   : temperature (eV)
C       Epar: parallel electric field (V/m)
C
C     Output:
C       Gamm: flux of primary runaway electrons (1/s)
C          
C----------------------------------------------------------------------|
      include 'src/inc/parameters_CODE_neural_network.inc'

        ! ----- Function arguments ------------------------------------|
      integer, intent(in) :: 
     >  nSpecies, Z(nSpecies), Z0(nSpecies)

      real(kind=dp), intent(in) ::
     >  Epar, nj(nSpecies), T

        ! ----- Local variables ---------------------------------------|
      real(kind=dp) ::
     >  input(8), ne, ne_tot

C----------------------------------------------------------------------|
C     Construct required quantities
C----------------------------------------------------------------------|
        ! ----- Free electron density ---------------------------------|
      ne   = sum(Z0(:)*nj(:))

        ! ----- Total electron density --------------------------------|
      ne_tot = sum(Z( :)*nj(:)) 

      input = (/    
     >  sum(Z0(:)**2*nj(:))/ne,
     >  sum((Z(:)**2-Z0(:)**2)*nj(:))/ne_tot,
     >  sum(nj(:)*Z0(:)/Z(:))/ne_tot,            
     >  sum(Z0(:)*Z(:)*nj(:))/ne_tot,
     >  log(ne), 
     >  ne/ne_tot,
     >  Epar/E_D(ne, T), 
     >  log(T/510998.946_dp)
     >  /)

        ! Normalize input
      input = (input - input_mean) / input_std

C----------------------------------------------------------------------|
C     Evaluate neural network
C----------------------------------------------------------------------|
      Gamma_D = 
     >  dot_product(W5, 
     >      tanh(matmul(W4, 
     >          tanh(matmul(W3, 
     >              tanh(matmul(W2, 
     >                  tanh(matmul(W1, input) + b1)) 
     >              + b2)) 
     >          + b3))
     >      + b4))
     >  + b5

      Gamma_D = 4._dp/(3*sqrt(pi))*nu_ee(ne, T)
     >  * exp(Gamma_D*output_std + output_mean)

      return
      end function Dreicer_growthrate_CODE_neural_network
C======================================================================|


      end module calculate_Dreicer_growthrate
C######################################################################|
