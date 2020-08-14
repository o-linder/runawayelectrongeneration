C######################################################################|
      module Coulomb_logarithms
C----------------------------------------------------------------------|
C     This module is used to calculate the Coulomb logarithm
C     associated with collisions of various particles.
C   
C     Functions:
C       - ln_Lambda_0(ne, T):       thermal <-> thermal electrons
C       - ln_Lambda_c(ne, T):       relativistic <-> relativistic 
C                                   electrons
C       - ln_Lambda_ee(ne, T, p):   thermal <-> non-thermal electrons
C       - ln_Lambda_ei(ne, T, p):   thermal ions <-> non-thermal
C                                   electrons
C
C     Function arguments:
C       - ne:   free electron density (m**-3)
C       - T:    electron temperature (eV)
C       - p:    normalised electron momentum
C
C----------------------------------------------------------------------|
      use double
      use physical_constants, only :
     >  m_ec

      implicit none

      interface ln_Lambda_0
          module procedure ln_Lambda_0_ra, ln_Lambda_0_rs
      end interface

      public ::
     >  ln_Lambda_0,
     >  ln_Lambda_c,
     >  ln_Lambda_ee,
     >  ln_Lambda_ei

      private ::
     >  ln_Lambda_0_ra, ln_Lambda_0_rs,
     >
     >  dp, k, m_ec

        ! ----- Parameters --------------------------------------------|
      real(kind=dp), parameter ::
     >  k   = 5._dp

      contains

C======================================================================|
C     Function ln_Lambda_0
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
      !----------------------------------------------------------------|
      function ln_Lambda_0_ra(ne, T)
      !----------------------------------------------------------------|
      real(kind=dp), dimension(:), intent(in) :: ne, T
      real(kind=dp), dimension(:), allocatable :: ln_Lambda_0_ra
      integer :: n

        ! Allocate array
      n = size(ne)
      allocate(ln_Lambda_0_ra(1:n))

        ! Calculate Coulomb logarithm
      ln_Lambda_0_ra(:) = 14.9_dp - .5_dp*log(1.e-20_dp*ne(:)) 
     >  + log(T(:)/1.e3_dp) 

      return
      end function ln_Lambda_0_ra
      !----------------------------------------------------------------|


      !----------------------------------------------------------------|
      real(kind=dp) function ln_Lambda_0_rs(ne, T)
      !----------------------------------------------------------------|
      real(kind=dp), intent(in) :: ne, T
      real(kind=dp), dimension(1) :: tmp

      tmp = ln_Lambda_0_ra((/ne/), (/T/))
      ln_Lambda_0_rs = tmp(1)

      return
      end function ln_Lambda_0_rs
      !----------------------------------------------------------------|
C======================================================================|


C======================================================================|
      real(kind=dp) function ln_Lambda_c(ne, T)
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
        ! ----- Function arguments ------------------------------------|
      real(kind=dp), intent(in) ::
     >  ne, T

C----------------------------------------------------------------------|
      ln_Lambda_c = 14.6_dp + .5_dp*log(T/(1.e-20_dp*ne)) 

      return
      end function ln_Lambda_c
C======================================================================|


C======================================================================|
      real(kind=dp) function ln_Lambda_ee(ne, T, p) 
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
        ! ----- Function arguments ------------------------------------|
      real(kind=dp), intent(in) ::
     >  ne, p, T

C----------------------------------------------------------------------|
      ln_Lambda_ee= ln_Lambda_0(ne, T) + 1._dp/k
     >  * log(1 + sqrt(m_ec*(sqrt(1+p**2)-1)/T)**k)

      return
      end function ln_Lambda_ee
C======================================================================|


C======================================================================|
      real(kind=dp) function ln_Lambda_ei(ne, T, p) 
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
        ! ----- Function arguments ------------------------------------|
      real(kind=dp), intent(in) ::
     >  ne, p, T

C----------------------------------------------------------------------|
      ln_Lambda_ei = ln_Lambda_0(ne, T) 
     >  + log(1 + (sqrt(2*m_ec/T)*p)**k)/k

      return
      end function ln_Lambda_ei
C======================================================================|


      end module Coulomb_logarithms
C######################################################################|
