C######################################################################|
      module collision_frequencies
C----------------------------------------------------------------------|
C     This module is used to calculate various collision frequencies.
C   
C     Functions:
C       - nu_ee(ne, T):             thermal electron-electron collisions
C       - nuBar_D(Z, Z0, nj, n...): normalised electron deflection 
C                                   frequency
C       - calc_nuBar_D01(Z, Z0...): first two terms of the expansion of
C                                   the normalised deflection frequency
C       - nuBar_S(Z, Z0, nj, n...): normalised electron slowing down 
C                                   frequency
C       - calc_nuBar_S01(Z, Z0...): first two terms of the expansion of
C                                   the normalised slowing down
C                                   frequency
C
C     Function and subroutine arguments:
C       - ne:   free electron density (m**-3)
C       - T:    electron temperature (eV)
C       - Z:    atomic number of each ion species
C       - Z0:   net charge of plasma for each ion species
C       - nj:   density of each ion species (m**-3)
C       - nSpe: number of ion species
C       - p:    normalised electron momentum
C
C----------------------------------------------------------------------|
      use double
      use physical_constants, only :
     >  alpha
      use Coulomb_logarithms

      implicit none

      public ::
     >  nu_ee,
     >  nuBar_D,
     >  calc_nuBar_D01,
     >  nuBar_S,
     >  calc_nuBar_S01

      private ::
     >  get_length_scale,
     >  get_Ihat,
     >
     >  alpha, dp, k, ln_Lambda_0, ln_Lambda_c, ln_Lambda_ee, 
     >  ln_Lambda_ei

        ! ----- Parameters --------------------------------------------|
      real(kind=dp), parameter ::
     >  k = 5.0_dp

      contains

C======================================================================|
      real(kind=dp) function nu_ee(ne, T)
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
C       nu_ee: thermal electron electron collision frequency
C
C----------------------------------------------------------------------|
        ! ----- Function arguments ------------------------------------|
      real(kind=dp), intent(in) ::
     >  ne, T

C----------------------------------------------------------------------|
      nu_ee = 3.863484401810658e-12_dp*ne*ln_Lambda_0(ne, T)/T**1.5_dp

      return
      end function nu_ee
C======================================================================|


C======================================================================|
      real(kind=dp) function nuBar_D(Z, Z0, nj, nSpecies, T, p) 
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
C
C----------------------------------------------------------------------|
        ! ----- Function arguments ------------------------------------|
      integer, intent(in) ::
     >  nSpecies, Z(nSpecies), Z0(nSpecies)

      real(kind=dp), intent(in) ::
     >  nj(nSpecies), p, T 

        ! ----- Local variables ---------------------------------------|
      integer ::
     >  i

      real(kind=dp) ::
     >  aBar, ne, paBar, Zeff

C----------------------------------------------------------------------|
C     Construct required quantities
C----------------------------------------------------------------------|
        ! ----- Free electron density ---------------------------------|
      ne = sum(Z0(:)*nj(:))    

        ! ----- Effective charge --------------------------------------|
      Zeff = sum(Z0(:)**2*nj(:))/ne

C----------------------------------------------------------------------|
C     Electron-ion contributions
C----------------------------------------------------------------------|
      nuBar_D = Zeff*ln_Lambda_ei(ne, T, p)

      do i=1,nSpecies
        aBar = get_length_scale(Z(i),Z0(i))

        if (aBar .gt. 0._dp) then
            paBar = (p*aBar)**(1.5_dp)

            nuBar_D = nuBar_D + nj(i)/ne/1.5_dp
     >          *(+ (Z(i)**2-Z0(i)**2)*log(paBar + 1) 
     >            - (Z(i)-Z0(i))**2*paBar/(paBar + 1))
        endif
      enddo

C----------------------------------------------------------------------|
C     Free electron-electron contributions
C----------------------------------------------------------------------|
      nuBar_D = (nuBar_D + ln_Lambda_ee(ne, T, p))/ln_Lambda_c(ne, T)

      return
      end function nuBar_D
C======================================================================|


C======================================================================|
      subroutine calc_nuBar_D01(Z, Z0, nj, nSpecies, T, nuBar_D0, 
     >  nuBar_D1)
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
C
C----------------------------------------------------------------------|
        ! ----- Subroutine arguments ----------------------------------|
      integer, intent(in) :: 
     >  nSpecies, Z(nSpecies), Z0(nSpecies)

      real(kind=dp), intent(in) ::
     >  nj(nSpecies), T

      real(kind=dp), intent(out) ::
     >  nuBar_D0, nuBar_D1

        ! ----- Local variables ---------------------------------------|
      integer ::
     >  i
      real(kind=dp) ::
     >  aBar, lnLambdac, ne, Zeff

C----------------------------------------------------------------------|
C     Construct required quantities
C----------------------------------------------------------------------|
        ! ----- Free electron density ---------------------------------|
      ne = sum(Z0(:)*nj(:))

        ! ----- Thermal Coulomb logarithm -----------------------------|
      lnLambdac = ln_Lambda_c(ne, T)

        ! ----- Effective charge --------------------------------------|
      Zeff = sum(Z0(:)**2*nj(:))/ne
      
C----------------------------------------------------------------------|
C     Construct expansion of generalized deflection frequency
C----------------------------------------------------------------------|
      nuBar_D0 = 1 + Zeff - sum((Z(:)-Z0(:))**2*nj(:))/1.5_dp/ne
     >  / lnLambdac

      do i=1,nSpecies
        aBar = get_length_scale(Z(i),Z0(i))

        if (aBar .gt. 0._dp) then
            nuBar_D0 = nuBar_D0 + nj(i)/ne/lnLambdac*
     >          (Z(i)**2-Z0(i)**2)*log(aBar)
        endif
      enddo

      nuBar_D1 = sum(Z(:)**2*nj(:))/ne/lnLambdac

      return
      end subroutine calc_nuBar_D01
C======================================================================|


C======================================================================|
      real(kind=dp) function nuBar_S(Z, Z0, nj, nSpecies, T, p) 
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
C
C----------------------------------------------------------------------|
        ! ----- Function arguments ------------------------------------|
      integer, intent(in) ::
     >  nSpecies, Z(nSpecies), Z0(nSpecies)

      real(kind=dp), intent(in) ::
     >  nj(nSpecies), p, T 

        ! ----- Local variables ---------------------------------------|
      integer ::
     >  i

      real(kind=dp) ::
     >  Ihat, ne, psg

C----------------------------------------------------------------------|
C     Construct required quantities
C----------------------------------------------------------------------|
        ! ----- Free electron density ---------------------------------|
      ne = sum(Z0(:)*nj(:))

        ! ----- Factor p*sqrt(gamma - 1)) -----------------------------|
      psg = p*sqrt(sqrt(p**2 + 1) - 1)

C----------------------------------------------------------------------|
C     Construct slowing-down frequency
C----------------------------------------------------------------------|
        ! ----- Add bound electron-electron contributions -------------|
      nuBar_S    = -sum(nj(:)*(Z(:)-Z0(:)))*p**2/(1 + p**2)

      do i=1,nSpecies
        Ihat    = get_Ihat(Z(i),Z0(i))

        if (Ihat .gt. 0._dp) 
     >      nuBar_S  = nuBar_S + nj(i)*(Z(i)-Z0(i))/k
     >          * log(1+(psg/Ihat)**k)
      enddo
      nuBar_S = nuBar_S/ne/ln_Lambda_c(ne, T)

        ! ----- Add free electron-electron contributions --------------|
      nuBar_S = nuBar_S + ln_Lambda_ee(ne, T, p)/ln_Lambda_c(ne, T)

      return
      end function nuBar_S
C======================================================================|


C======================================================================|
      subroutine calc_nuBar_S01(Z, Z0, nj, nSpecies, T, nuBar_S0, 
     >  nuBar_S1)
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
C
C----------------------------------------------------------------------|
        ! ----- Subroutine arguments ----------------------------------|
      integer, intent(in) :: 
     >  nSpecies, Z(nSpecies), Z0(nSpecies)

      real(kind=dp), intent(in) ::
     >  nj(nSpecies), T

      real(kind=dp), intent(out) ::
     >  nuBar_S0, nuBar_S1

        ! ----- Local variables ---------------------------------------|
      integer ::
     >  i

      real(kind=dp) ::
     >  ne

C----------------------------------------------------------------------|
C     Construct required quantities
C----------------------------------------------------------------------|
        ! ----- Free electron density ---------------------------------|
      ne = sum(Z0(:)*nj(:))

C----------------------------------------------------------------------|
C     Construct expansion of slowing-down frequency
C----------------------------------------------------------------------|
      nuBar_S0 = sum((Z(:)-Z0(:))*nj(:))
      do i=1,nSpecies
        nuBar_S0 = nuBar_S0 + (Z(i)-Z0(i))*nj(i)
     >      * log(get_Ihat(Z(i),Z0(i)))
      enddo
      nuBar_S0 = 1 - nuBar_S0/ne/ln_Lambda_c(ne, T)

      nuBar_S1 = (1.5_dp*sum(Z(:)*nj(:))/ne - 1)/ln_Lambda_c(ne, T)

      return
      end subroutine calc_nuBar_S01
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
C
C----------------------------------------------------------------------|
      include 'src/inc/parameters_get_length_scale.inc'

        ! ----- Function arguments ------------------------------------|
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

      return
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
C
C----------------------------------------------------------------------|
      include 'src/inc/parameters_get_Ihat.inc'

        ! ----- Function arguments ------------------------------------|
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

      return
      end function get_Ihat
C======================================================================|


      end module collision_frequencies
C######################################################################|
