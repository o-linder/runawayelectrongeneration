C######################################################################|
      module runawayElectronGeneration

      use double

      implicit none

      public  :: 
            ! RE growth rates
     >  calc_Gamma_av,
     >  calc_Gamma_av_RP,
     >  calc_Gamma_D,
     >  calc_Gamma_D_nn,
     >
            ! Electric fields
     >  calc_EceffOverEc,
     >  calc_Ec,
     >  calc_ED,
     >
            ! Critical momentum
     >  calc_pstar,
     >
            ! Collision frequencies
     >  calc_nuBarD,
     >  calc_nuBarD01,
     >  calc_nuBarS,
     >  calc_nuBarS01,
     >  calc_nuee,
     >
            ! Coulomb logarithms
     >  calc_lnLambda0,
     >  calc_lnLambdac,
     >  calc_lnLambdaee,
     >  calc_lnLambdaei

      private ::
     >  get_length_scale,
     >  get_Ihat,
     >
     >  alpha, dp, iter_acc, iter_max, pi

        ! ----- Parameters --------------------------------------------|
      integer, parameter ::
     >  iter_max    = 100               ! maximum number of iterations

      real(kind=dp), parameter ::
     >  alpha       = 1/137.036_dp,
     >  iter_acc    = 1.0e-8_dp,        ! target accuracy of iterations
     >  pi          = 3.1415926535897931_dp

      contains

C======================================================================|
      real(kind=dp) function calc_Gamma_av_RP(Z, Z0, nj, nSpecies, T, 
     >  Epar, eps) result (Gamma_av)
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
C----------------------------------------------------------------------|
        ! ----- Function variables ------------------------------------|
      integer, intent(in) :: 
     >  nSpecies, Z(nSpecies), Z0(nSpecies)
      real(kind=dp), intent(in) ::
     >  Epar, eps, nj(nSpecies), T

        ! ----- Local variables ---------------------------------------|
      integer ::
     >  i
      real(kind=dp) ::
     >  Ec, gam, lnLambdac, ne_free, nuee, Zeff

C----------------------------------------------------------------------|
C     Construct required quantities
C----------------------------------------------------------------------|
        ! ----- Free electron density ---------------------------------|
      ne_free   = sum(Z0(:)*nj(:))

        ! ----- Critical electric field -------------------------------|
      Ec        = calc_Ec(ne_free, T)

        ! ----- Relativistic Coulomb logarithm ------------------------|
      lnLambdac = calc_lnLambdac(ne_free, T) 

        ! ----- Thermal electron electron collision frequency ---------|
      nuee      = calc_nuee(ne_free, T)
     >              *calc_lnLambda0(ne_free,T)/lnLambdac
     >              *7.7431021134304697e-09_dp*T**1.5_dp

        ! ----- Effective charge --------------------------------------|
      Zeff      = sum(Z0(:)**2*nj(:))/ne_free

C----------------------------------------------------------------------|
C     Calculate growth rate
C----------------------------------------------------------------------|
      gam       = 1._dp/(1 + 1.46_dp*sqrt(eps) + 1.72_dp*eps)

      Gamma_av  = nuee*sqrt(pi*gam/3._dp/(Zeff+5))*(Epar/Ec-1)/lnLambdac
     >          / sqrt(1 - Ec/Epar + 4*pi*(Zeff+1)**2/3._dp/gam/(Zeff+5)
     >          / (Epar**2/Ec**2 + 4._dp/gam**2-1))

      end function calc_Gamma_av_RP
C======================================================================|

C======================================================================|
      real(kind=dp) function calc_Gamma_av(Z, Z0, nj, nSpecies, T, 
     >  Epar, B) result (Gamma_av) 
C----------------------------------------------------------------------|
C     The flux of secondary electrons to relativistic momenta is 
C     calculated according to Eq. (14) in
C     L. Hesslow et al., Nucl. Fusion 59, 084004 (2019).
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
C----------------------------------------------------------------------|
        ! ----- Function variables ------------------------------------|
      integer, intent(in) :: 
     >  nSpecies, Z(nSpecies), Z0(nSpecies)
      real(kind=dp), intent(in) ::
     >  B, Epar, nj(nSpecies), T

        ! ----- Local variables ---------------------------------------|
      integer ::
     >  i
      real(kind=dp) ::
     >  Ec, lnLambdac, ne_free, ne_tot, pstar

C----------------------------------------------------------------------|
C     Construct required quantities
C----------------------------------------------------------------------|
        ! ----- Free electron density ---------------------------------|
      ne_free   = sum(Z0(:)*nj(:))

        ! ----- Critical electric field -------------------------------|
      Ec        = calc_Ec(ne_free, T)

        ! ----- Relativistic Coulomb logarithm ------------------------|
      lnLambdac = calc_lnLambdac(ne_free, T) 

        ! ----- Total electron density --------------------------------|
      ne_tot    = sum(Z( :)*nj(:)) 

        ! ----- Effective critical momentum for runaway ---------------| 
      pstar     = calc_pstar(Z, Z0, nj, nSpecies, T, Epar, B)

C----------------------------------------------------------------------|
C     Calculate growth rate
C----------------------------------------------------------------------|
      Gamma_av  = 5.8667920978923496e+02_dp*ne_tot*Ec/ne_free/lnLambdac
     >          * (Epar/Ec - calc_EceffOverEc(Z, Z0, nj, nSpecies, T,B))
     >          / sqrt(4 + calc_nuBarD(Z, Z0, nj, nSpecies, T, pstar)
     >                  *  calc_nuBarS(Z, Z0, nj, nSpecies, T, pstar)) 

      end function calc_Gamma_av
C======================================================================|

C======================================================================|
      real(kind=dp) function calc_Gamma_D(Z, Z0, nj, nSpecies, T, 
     >  Epar) result (Gamma_D) 
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
C----------------------------------------------------------------------|
        ! ----- Parameters --------------------------------------------|
      real(kind=dp), parameter ::
     >  CR = 1.0_dp

        ! ----- Function variables ------------------------------------|
      integer, intent(in) :: 
     >  nSpecies, Z(nSpecies), Z0(nSpecies)
      real(kind=dp), intent(in) ::
     >  Epar, nj(nSpecies), T

        ! ----- Local variables ---------------------------------------|
      integer ::
     >  i
      real(kind=dp) ::
     >  Ec, ED, gam, h, lambda, ne_free, Zeff

C----------------------------------------------------------------------|
C     Construct required quantities
C----------------------------------------------------------------------|
        ! ----- Free electron density ---------------------------------|
      ne_free   = sum(Z0(:)*nj(:))

        ! ----- Critical electric field -------------------------------|
      Ec        = Epar/calc_Ec(ne_free, T)

        ! ----- Dreicer electric field --------------------------------|
      ED        = Epar/calc_ED(ne_free, T)

        ! ----- Effective charge --------------------------------------|
      Zeff      = sum(Z0(:)**2*nj(:))/ne_free

C----------------------------------------------------------------------|
C     Calculate growth rate
C----------------------------------------------------------------------|
      Gamma_D   = 0._dp
      if (Ec .le. 1._dp) return

      h         = (Ec - (Zeff-7)/(Zeff+1) + 2*sqrt(Ec/(Ec-1))*(Ec-2))
     >              *(Zeff+1)/(16*(Ec-1))

      lambda    = 8*Ec**2*(1 - 1._dp/(2*Ec) - sqrt(1 - 1._dp/Ec))

      gam       = Ec/2*sqrt((Zeff+1)/(Ec-1))
     >              *(3.141592653589793_dp/2 - asin(1-2/Ec))

      Gamma_D   = CR*calc_nuee(ne_free, T)*ED**(-h)
     >              *exp(-lambda/(4*ED) - gam/sqrt(ED))

      end function calc_Gamma_D
C======================================================================|

C======================================================================|
      real(kind=dp) function calc_Gamma_D_nn(Z, Z0, nj, nSpecies, T, 
     >  Epar) result (Gamma_D)
C----------------------------------------------------------------------|
C     The flux of primary electrons to relativistic momenta is 
C     calculated according to L. Hesslow et al., To be submitted to
C     J. Plasma Phys.
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
C----------------------------------------------------------------------|
      include 'parameters_calc_Gamma_D_nn.inc'

        ! ----- Function variables ------------------------------------|
      integer, intent(in) :: 
     >  nSpecies, Z(nSpecies), Z0(nSpecies)
      real(kind=dp), intent(in) ::
     >  Epar, nj(nSpecies), T

        ! ----- Local variables ---------------------------------------|
      real(kind=dp) ::
     >  input(8), ne_free, ne_tot

C----------------------------------------------------------------------|
C     Construct required quantities
C----------------------------------------------------------------------|
        ! ----- Free electron density ---------------------------------|
      ne_free   = sum(Z0(:)*nj(:))

        ! ----- Total electron density --------------------------------|
      ne_tot    = sum(Z( :)*nj(:)) 

      input     = (/    
     >      sum(Z0(:)**2*nj(:))/ne_free,
     >      sum((Z(:)**2-Z0(:)**2)*nj(:))/ne_tot,
     >      sum(nj(:)*Z0(:)/Z(:))/ne_tot,            
     >      sum(Z0(:)*Z(:)*nj(:))/ne_tot,
     >      log(ne_free), 
     >      ne_free/ne_tot,
     >      Epar/calc_ED(ne_free, T), 
     >      log(T/510998.946_dp)
     >      /)

        ! Normalize input
      input = (input - input_mean) / input_std

C----------------------------------------------------------------------|
C     Evaluate neural network
C----------------------------------------------------------------------|
      Gamma_D   = 
     >      dot_product(W5, 
     >          tanh(matmul(W4, 
     >              tanh(matmul(W3, 
     >                  tanh(matmul(W2, 
     >                      tanh(matmul(W1, input) + b1)) 
     >                  + b2)) 
     >              + b3))
     >          + b4))
     >      + b5

      Gamma_D   = 4._dp/(3._dp*sqrt(pi))*calc_nuee(ne_free, T)
     >              * exp(Gamma_D*output_std + output_mean)

      end function calc_Gamma_D_nn
C======================================================================|

C======================================================================|
      real(kind=dp) function calc_EceffOverEc(Z, Z0, nj, nSpecies,
     >  T, B) result(EceffOverEc)
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
C----------------------------------------------------------------------|
        ! ----- Function variables ------------------------------------|
      integer, intent(in) ::
     >  nSpecies, Z(nSpecies), Z0(nSpecies)
      real(kind=dp), intent(in) ::
     >  B, nj(nSpecies), T

        ! ----- Local variables ---------------------------------------|
      integer ::
     >  i
      real(kind=dp) ::
     >  e0, e1, e2, Eold, ne_free, 
     >  nuBarD0, nuBarD1, nuBarS0, nuBarS1, 
     >  tauRadInv

C----------------------------------------------------------------------|
C     Construct required quantities
C----------------------------------------------------------------------|
        ! ----- Free electron density ---------------------------------|
      ne_free   = sum(Z0(:)*nj(:))
      
        ! ----- Generalised deflection frequency ----------------------|
      call calc_nuBarD01(Z, Z0, nj, nSpecies, T, nuBarD0, nuBarD1)

        ! ----- Slowing down frequency --------------------------------|
      call calc_nuBarS01(Z, Z0, nj, nSpecies, T, nuBarS0, nuBarS1)
  
        ! ----- Synchrotron radiation term ----------------------------|
      tauRadInv = B**2/(15.44_dp*calc_lnLambdac(ne_free, T))
     >              /(1.e-20_dp*ne_free)

C----------------------------------------------------------------------|
C     Calculate normalised Eceff iteratively
C----------------------------------------------------------------------|
        ! ----- Coefficients ------------------------------------------|
      e0    = nuBarS0 + nuBarS1 * (1+nuBarD1/nuBarD0)
     >          * log(nuBarD0/(2*nuBarS1))
      e1    = (.7_dp+.4_dp*log(nuBarD0/(2*nuBarS1)))
     >          *alpha*nuBarD0*nuBarD1 + nuBarS1**2
      e2    = 2*tauRadInv*nuBarD0**2

        ! ----- Iterative calculation ---------------------------------|
      EceffOverEc   = calc_lnLambdac(ne_free, T)
     >                  /calc_lnLambda0(ne_free, T)
     >                  * sum(Z(:)*nj(:))/ne_free

      do i=1,iter_max
        Eold        = EceffOverEc
        EceffOverEc = e0 + sqrt(e1 + e2/Eold)
        if (abs(EceffOverEc-Eold) .le. iter_acc) 
     >      exit 
      enddo

      end function calc_EceffOverEc
C======================================================================|

C======================================================================|
      real(kind=dp) function calc_Ec(ne, T) result(Ec)
C----------------------------------------------------------------------|
C     Calculates the critical electric field Ec:
C       Ec = me*c/e/tauc = ne*e**3*ln Lambda_c/(4*pi*eps0**2*me*c**2).
C
C     Input:
C       ne  : free electron density (m**-3)
C       T   : electron temperature (eV)
C
C     Output:
C       Ec  : critical electric field
C----------------------------------------------------------------------|
        ! ----- Function variables ------------------------------------|
      real(kind=dp), intent(in) ::
     >  ne, T

C----------------------------------------------------------------------|
      Ec    = 5.09909908817339515e-23_dp*ne*calc_lnLambdac(ne, T)

      end function calc_Ec
C======================================================================|

C======================================================================|
      real(kind=dp) function calc_ED(ne, T) result(ED)
C----------------------------------------------------------------------|
C     Calculates the Dreicer electric field ED:
C       ED = ne*e**3*ln Lambda/(4*pi*eps0**2*T).
C
C     Input:
C       ne  : free electron density (m**-3)
C       T   : electron temperature (eV)
C
C     Output:
C       ED  : Dreicer electric field
C----------------------------------------------------------------------|
        ! ----- Function variables ------------------------------------|
      real(kind=dp), intent(in) ::
     >  ne, T

C----------------------------------------------------------------------|
      ED    = 2.6056342609758366e-17_dp*ne*calc_lnLambda0(ne, T)/T

      end function calc_ED
C======================================================================|

C======================================================================|
      real(kind=dp) function calc_pstar(Z, Z0, nj, nSpecies,
     >  T, Epar, B) result(pstar)
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
C----------------------------------------------------------------------|
        ! ----- Function variables ------------------------------------|
      integer, intent(in) ::
     >  nSpecies, Z(nSpecies), Z0(nSpecies)
      real(kind=dp), intent(in) ::
     >  B, nj(nSpecies), T, Epar 

        ! ----- Local variables ---------------------------------------|
      integer ::
     >  i
      real(kind=dp) ::
     >  EEc, ne_free, pold

C----------------------------------------------------------------------|
C     Construct required quantities
C----------------------------------------------------------------------|
        ! ----- Free electron density ---------------------------------|
      ne_free   = sum(Z0(:)*nj(:))

        ! ----- Normalised electric field -----------------------------|
      EEc       = max(  Epar/calc_Ec(ne_free, T),
     >                  calc_EceffOverEc(Z, Z0, nj, nSpecies, T, B))

C----------------------------------------------------------------------|
C     Obtain critical momentum iteratively
C----------------------------------------------------------------------|
      pstar     = 1._dp

      do i=1,iter_max
        pold    = pstar
        pstar   = sqrt(1._dp/EEc*sqrt(
     >        calc_nuBarS(Z(:), Z0(:), nj(:), nSpecies, T, pstar)
     >      * calc_nuBarD(Z(:), Z0(:), nj(:), nSpecies, T, pstar)))
        if (abs(pstar-pold) .le. iter_acc)
     >      exit 
      enddo

      end function calc_pstar
C======================================================================|

C======================================================================|
      real(kind=dp) function calc_nuBarD(Z, Z0, nj, nSpecies, T, p) 
     >  result(nuBarD)
C----------------------------------------------------------------------|
C     Calculates the electron deflection frequency, normalised to the 
C     relativistic electron collision time according to Eq. (2.22) in 
C     L. Hesslow et al., J. Plasma Phys. 84, 905840605 (2018) and Eq.
C     (5) in L. Hesslow et al., Plasma Phys. Control. Fusion 60, 074010
C     (2018).
C
C     Input:
C       Z   : atomic number of each ion species
C       Z0  : net charge of plasma for each ion species
C       nj  : density of each ion species (m**-3)
C       nSpe: number of ion species
C           Z, Z0 and nj must be vectors of length nSpecies
C       T   : temperature (eV)
C       p   : normalised electron momentum
C
C     Output:
C       nuBD: electron deflection frequency
C----------------------------------------------------------------------|
        ! ----- Function variables ------------------------------------|
      integer, intent(in) ::
     >  nSpecies, Z(nSpecies), Z0(nSpecies)
      real(kind=dp), intent(in) ::
     >  nj(nSpecies), p, T 

        ! ----- Local variables ---------------------------------------|
      integer ::
     >  i
      real(kind=dp) ::
     >  aBar, lnLambdac, ne_free, paBar, Zeff

C----------------------------------------------------------------------|
C     Construct required quantities
C----------------------------------------------------------------------|
        ! ----- Free electron density ---------------------------------|
      ne_free   = sum(Z0(:)*nj(:))    

        ! ----- Relativistic Coulomb logarithm ------------------------|
      lnLambdac = calc_lnLambdac(ne_free, T)

        ! ----- Effective charge --------------------------------------|
      Zeff      = sum(Z0(:)**2*nj(:))/ne_free

C----------------------------------------------------------------------|
C     Electron-ion contributions
C----------------------------------------------------------------------|
      nuBarD    = Zeff*calc_lnLambdaei(ne_free, T, p)

      do i=1,nSpecies
        aBar    = get_length_scale(Z(i),Z0(i))
        if (aBar .gt. 0._dp) then
            paBar   = (p*aBar)**(1.5_dp)
            nuBarD  = nuBarD + nj(i)/ne_free/1.5_dp*(
     >                    (Z(i)**2-Z0(i)**2)*log(paBar + 1) 
     >                  - (Z(i)-Z0(i))**2*paBar/(paBar + 1))
        endif
      enddo

C----------------------------------------------------------------------|
C     Free electron-electron contributions
C----------------------------------------------------------------------|
      nuBarD    = (nuBarD + calc_lnLambdaee(ne_free, T, p))/lnLambdac

      end function calc_nuBarD
C======================================================================|

C======================================================================|
      subroutine calc_nuBarD01(Z, Z0, nj, nSpecies, T, nuBarD0, nuBarD1) 
C----------------------------------------------------------------------|
C     Calculates the first two terms of the expansion of the generalized
C     electron deflection frequency in ln p as in L. Hesslow et al., 
C     Plasma Phys. Control. Fusion 60, 074010 (2018) in Eqs. (6)-(8).
C
C       Input:
C           Z   : atomic number of each ion species
C           Z0  : net charge of plasma for each ion species
C           nj  : density of each ion species (m**-3)
C           nSpe: number of ion species
C               Z, Z0 and nj must be vectors of length nSpecies
C           T   : temperature (eV)
C
C       Output:
C           nBD0: 0th order electron deflection frequency
C           nBD1: 1st order electron deflection frequency
C----------------------------------------------------------------------|
        ! ----- Subroutine variables ----------------------------------|
      integer, intent(in) :: 
     >  nSpecies, Z(nSpecies), Z0(nSpecies)
      real(kind=dp), intent(in) ::
     >  nj(nSpecies), T

      real(kind=dp), intent(out) ::
     >  nuBarD0, nuBarD1

        ! ----- Local variables ---------------------------------------|
      integer ::
     >  i
      real(kind=dp) ::
     >  aBar, lnLambdac, ne_free, Zeff

C----------------------------------------------------------------------|
C     Construct required quantities
C----------------------------------------------------------------------|
        ! ----- Free electron density ---------------------------------|
      ne_free   = sum(Z0(:)*nj(:))

        ! ----- Thermal Coulomb logarithm -----------------------------|
      lnLambdac = calc_lnLambdac(ne_free, T)

        ! ----- Effective charge --------------------------------------|
      Zeff      = sum(Z0(:)**2*nj(:))/ne_free
      
C----------------------------------------------------------------------|
C     Construct expansion of generalized deflection frequency
C----------------------------------------------------------------------|
      nuBarD0   = 1 + Zeff - sum((Z(:)-Z0(:))**2*nj(:))
     >              /1.5_dp/ne_free/lnLambdac
      do i=1,nSpecies
        aBar    = get_length_scale(Z(i),Z0(i))
        if (aBar .gt. 0._dp) 
     >      nuBarD0 = nuBarD0 + nj(i)/ne_free/lnLambdac*
     >                  (Z(i)**2-Z0(i)**2)*log(aBar)
      enddo

      nuBarD1   = sum(Z(:)**2*nj(:))/ne_free/lnLambdac

      end subroutine calc_nuBarD01
C======================================================================|

C======================================================================|
      real(kind=dp) function calc_nuBarS(Z, Z0, nj, nSpecies, T, p) 
     >  result(nuBarS)
C----------------------------------------------------------------------|
C     Calculates the electron slowing down frequency, normalised to the 
C     relativistic electron collision time according to Eq. (2.31) in 
C     L. Hesslow et al., J. Plasma Phys. 84, 905840605 (2018).
C
C     Input:
C       Z   : atomic number of each ion species
C       Z0  : net charge of plasma for each ion species
C       nj  : density of each ion species (m**-3)
C       nSpe: number of ion species
C           Z, Z0 and nj must be vectors of length nSpecies
C       T   : temperature (eV)
C       p   : normalised electron momentum
C
C     Output:
C       nuBS: electron slowing down frequency
C----------------------------------------------------------------------|
        ! ----- Parameters --------------------------------------------|
      real(kind=dp), parameter ::
     >  k   = 5.0_dp

        ! ----- Function variables ------------------------------------|
      integer, intent(in) ::
     >  nSpecies, Z(nSpecies), Z0(nSpecies)
      real(kind=dp), intent(in) ::
     >  nj(nSpecies), p, T 

        ! ----- Local variables ---------------------------------------|
      integer ::
     >  i
      real(kind=dp) ::
     >  Ihat, lnLambdac, ne_free, psg

C----------------------------------------------------------------------|
C     Construct required quantities
C----------------------------------------------------------------------|
        ! ----- Free electron density ---------------------------------|
      ne_free   = sum(Z0(:)*nj(:))

        ! ----- Relativistic Coulomb logarithm ------------------------|
      lnLambdac = calc_lnLambdac(ne_free, T)

        ! ----- Factor p*sqrt(gamma - 1)) -----------------------------|
      psg       = p*sqrt(sqrt(p**2 + 1) - 1)

C----------------------------------------------------------------------|
C     Construct slowing-down frequency
C----------------------------------------------------------------------|
        ! ----- Add bound electron-electron contributions -------------|
      nuBarS    = -sum(nj(:)*(Z(:)-Z0(:)))*p**2/(1 + p**2)

      do i=1,nSpecies
        Ihat    = get_Ihat(Z(i),Z0(i))
        if (Ihat .gt. 0._dp) 
     >      nuBarS  = nuBarS + nj(i)*(Z(i)-Z0(i))/k*log(1+(psg/Ihat)**k)
      enddo
      nuBarS        = nuBarS/ne_free/lnLambdac

        ! ----- Add free electron-electron contributions --------------|
      nuBarS       = nuBarS + calc_lnLambdaee(ne_free, T, p)/lnLambdac

      end function calc_nuBarS
C======================================================================|

C======================================================================|
      subroutine calc_nuBarS01(Z, Z0, nj, nSpecies, T, nuBarS0, nuBarS1)
C----------------------------------------------------------------------|
C     Calculates the first two terms of the expansion of the electron 
C     slowing-down frequency, normalized to the relativistic electron
C     collision time,  in ln p as in L. Hesslow et al., Plasma Phys. 
C     Control. Fusion 60, 074010 (2018) in Eqs. (10)-(12).
C
C       Input:
C           Z   : atomic number of each ion species
C           Z0  : net charge of plasma for each ion species
C           nj  : density of each ion species (m**-3)
C           nSpe: number of ion species
C               Z, Z0 and nj must be vectors of length nSpecies
C           T   : temperature (eV)
C
C       Output:
C           nBS0: 0th order electron slowing-down frequency
C           nBS1: 1st order electron slowing-down frequency
C----------------------------------------------------------------------|
        ! ----- Subroutine variables ----------------------------------|
      integer, intent(in) :: 
     >  nSpecies, Z(nSpecies), Z0(nSpecies)
      real(kind=dp), intent(in) ::
     >  nj(nSpecies), T

      real(kind=dp), intent(out) ::
     >  nuBarS0, nuBarS1

        ! ----- Local variables ---------------------------------------|
      integer ::
     >  i
      real(kind=dp) ::
     >  lnLambdac, ne_free

C----------------------------------------------------------------------|
C     Construct required quantities
C----------------------------------------------------------------------|
        ! ----- Free electron density ---------------------------------|
      ne_free   = sum(Z0(:)*nj(:))

        ! ----- Relativistic Coulomb logarithm ------------------------|
      lnLambdac = calc_lnLambdac(ne_free, T)

C----------------------------------------------------------------------|
C     Construct expansion of slowing-down frequency
C----------------------------------------------------------------------|
      nuBarS0   = sum((Z(:)-Z0(:))*nj(:))
      do i=1,nSpecies
        nuBarS0 = nuBarS0 +(Z(i)-Z0(i))*nj(i)*log(get_Ihat(Z(i),Z0(i)))
      enddo
      nuBarS0   = 1 - nuBarS0/ne_free/lnLambdac

      nuBarS1   = (1.5_dp*sum(Z(:)*nj(:))/ne_free - 1)/lnLambdac

      end subroutine calc_nuBarS01
C======================================================================|

C======================================================================|
      real(kind=dp) function calc_nuee(ne, T) result(nuee)
C----------------------------------------------------------------------|
C     Calculates the thermal electron-electron collision frequency
C     according to P. Helander et al., Plasma Phys. Control. Fusion 44,
C     B247 (2002).
C
C     Input:
C       ne  : free electron density (m**-3)
C       T   : electron temperature (eV)
C
C     Output:
C       nuee: thermal electron electron collision frequency
C----------------------------------------------------------------------|
        ! ----- Function variables ------------------------------------|
      real(kind=dp), intent(in) ::
     >  ne, T

C----------------------------------------------------------------------|
      nuee = 3.863484401810658e-12_dp*ne*calc_lnLambda0(ne, T)/T**1.5_dp

      end function calc_nuee
C======================================================================|

C======================================================================|
      real(kind=dp) function calc_lnLambda0(ne, T) result(lnLambda0)
C----------------------------------------------------------------------|
C     Calculates the thermal electron Coulomb logarithm according to 
C     Eq. (2.7) in L. Hesslow et al., J. Plasma Phys. 84, 905840605 
C     (2018).
C
C     Input:
C       ne  : free electron density (m**-3)
C       T   : electron temperature (eV)
C
C     Output:
C       lnL0: thermal electron Coulomb logarithm
C----------------------------------------------------------------------|
        ! ----- Function variables ------------------------------------|
      real(kind=dp), intent(in) ::
     >  ne, T

C----------------------------------------------------------------------|
      lnLambda0 = 14.9_dp - 0.5_dp*log(1.0e-20_dp*ne) + log(T/1.0e3_dp) 

      end function calc_lnLambda0
C======================================================================|

C======================================================================|
      real(kind=dp) function calc_lnLambdac(ne, T) result(lnLambdac)
C----------------------------------------------------------------------|
C     Calculates the relativistic electron Coulomb logarithm according 
C     to Eq. (2.9) in L. Hesslow et al., J. Plasma Phys. 84, 905840605 
C     (2018).
C
C     Input:
C       ne  : free electron density (m**-3)
C       T   : electron temperature (eV)
C
C     Output:
C       lnLc: relativistic electron Coulomb logarithm
C----------------------------------------------------------------------|
        ! ----- Function variables ------------------------------------|
      real(kind=dp), intent(in) ::
     >  ne, T

C----------------------------------------------------------------------|
      lnLambdac = 14.6_dp + 0.5_dp*log(T/(1.0e-20_dp*ne)) 

      end function calc_lnLambdac
C======================================================================|

C======================================================================|
      real(kind=dp) function calc_lnLambdaee(ne, T, p) 
     >  result(lnLambdaee)
C----------------------------------------------------------------------|
C     Calculates the Coulomb logarithm for collisions between
C     non-thermal and thermal electrons according to Eq. (2.10) in 
C     L. Hesslow et al., J. Plasma Phys. 84, 905840605 (2018).
C
C     Input:
C       ne  : free electron density (m**-3)
C       T   : electron temperature (eV)
C       p   : normalised electron momentum
C
C     Output:
C       lnLe: electron-electron Coulomb logarithm
C----------------------------------------------------------------------|
        ! ----- Parameters --------------------------------------------|
      real(kind=dp), parameter ::
     >  mc2 = 0.5109989461e6_dp,
     >  k   = 5._dp

        ! ----- Function variables ------------------------------------|
      real(kind=dp), intent(in) ::
     >  ne, p, T

C----------------------------------------------------------------------|
      lnLambdaee= calc_lnLambda0(ne, T) + 1._dp/k
     >      * log(1 + sqrt(mc2*(sqrt(1+p**2)-1)/T)**k)

      end function calc_lnLambdaee
C======================================================================|

C======================================================================|
      real(kind=dp) function calc_lnLambdaei(ne, T, p) 
     >  result(lnLambdaei)
C----------------------------------------------------------------------|
C     Calculates the Coulomb logarithm for collisions between
C     non-thermal electrons and ions according to Eq. (2.10) in 
C     L. Hesslow et al., J. Plasma Phys. 84, 905840605 (2018).
C
C     Input:
C       ne  : free electron density (m**-3)
C       T   : electron temperature (eV)
C       p   : normalised electron momentum
C
C     Output:
C       lnLi: electron-ion Coulomb logarithm
C----------------------------------------------------------------------|
        ! ----- Parameters --------------------------------------------|
      real(kind=dp), parameter ::
     >  mc2 = 0.5109989461e6_dp,
     >  k   = 5._dp

        ! ----- Function variables ------------------------------------|
      real(kind=dp), intent(in) ::
     >  ne, p, T

C----------------------------------------------------------------------|
      lnLambdaei    = calc_lnLambda0(ne, T) 
     >                  + log(1 + (sqrt(2*mc2/T)*p)**k)/k

      end function calc_lnLambdaei
C======================================================================|

C======================================================================|
      real(kind=dp) function get_length_scale(Z, Z0) result(aBar)
C----------------------------------------------------------------------|
C     Retrieves the effective length scale aBar = 2*a/alpha.
C     
C     For details, see L. Hesslow et al., J. Plasma Phys. 84, 905840605
C     (2018), particularly Table 1 and Appendix B.
C
C     Input:
C       Z   : atomic number of ion species
C       Z0  : net charge of ion species
C
C     Output:
C       a   : effective length scale
C----------------------------------------------------------------------|
      include 'parameters_get_length_scale.inc'

        ! ----- Function variables ------------------------------------|
      integer, intent(in) ::
     >  Z, Z0

C----------------------------------------------------------------------|
C     Effective length scale from the Kirillov model
C----------------------------------------------------------------------|
      aBar  = 0.7616184731724444_dp*(Z-Z0)**(2._dp/3._dp)/Z

C----------------------------------------------------------------------|
C     Effective length scale from density functional theory
C----------------------------------------------------------------------|
      if     (Z .eq.  2) then
        if (Z0 .lt. size(vHe))  aBar = vHe(Z0+1)
      elseif (Z .eq.  4) then
        if (Z0 .lt. size(vBe))  aBar = vBe(Z0+1)
      elseif (Z .eq.  6) then
        if (Z0 .lt. size(vC ))  aBar = vC( Z0+1)
      elseif (Z .eq.  7) then
        if (Z0 .lt. size(vN ))  aBar = vN( Z0+1)
      elseif (Z .eq. 10) then
        if (Z0 .lt. size(vNe))  aBar = vNe(Z0+1)
      elseif (Z .eq. 18) then
        if (Z0 .lt. size(vAr))  aBar = vAr(Z0+1)
      elseif (Z .eq. 54) then
        if (Z0 .lt. size(vXe))  aBar = vXe(Z0+1)
      endif

C----------------------------------------------------------------------|
C     Apply normalisation
C----------------------------------------------------------------------|
      aBar  = 2*aBar/alpha

      end function get_length_scale
C======================================================================|

C======================================================================|
      real(kind=dp) function get_Ihat(Z, Z0) result(Ihat) 
C----------------------------------------------------------------------|
C     Retrieves the normalized ion mean excitation energy I 
C     (normalized to electron rest mass).
C     
C     For details, see L. Hesslow et al., Plasma Phys. Control. Fusion 
C     60, 074010 (2018).
C
C     Input:
C       Z   : atomic number of ion species
C       Z0  : net charge of ion species
C
C     Output:
C       Ihat: normalised ion mean excitation energy
C----------------------------------------------------------------------|
      include 'parameters_get_Ihat.inc'

        ! ----- Function variables ------------------------------------|
      integer, intent(in) ::
     >  Z, Z0

C----------------------------------------------------------------------|
      Ihat   = 1._dp
      if (Z0 .lt. 0) return

      if     (Z .eq.  2) then
        if (Z0 .lt. size(vHe))  Ihat = vHe(Z0+1)
      elseif (Z .eq.  3) then
        if (Z0 .lt. size(vLi))  Ihat = vLi(Z0+1)
      elseif (Z .eq.  6) then
        if (Z0 .lt. size(vC ))  Ihat = vC (Z0+1)
      elseif (Z .eq. 10) then
        if (Z0 .lt. size(vNe))  Ihat = vNe(Z0+1)
      elseif (Z .eq. 18) then
        if (Z0 .lt. size(vAr))  Ihat = vAr(Z0+1)
      elseif (Z .eq. 36) then
        if (Z0 .lt. size(vKr))  Ihat = vKr(Z0+1)
      elseif (Z .eq. 54) then
        if (Z0 .lt. size(vXe))  Ihat = vXe(Z0+1)
      endif

      end function get_Ihat
C======================================================================|

      end module runawayElectronGeneration
C######################################################################|
