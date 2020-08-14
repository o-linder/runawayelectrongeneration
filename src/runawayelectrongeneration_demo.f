C======================================================================|
      program runawayElectronGeneration_demo

      use double
      use runawayelectrongeneration
      use calculate_hot_tail_population

      implicit none

      real(kind=dp) ::
     >  ne, p, Te, E, tmp

      real(kind=dp), dimension(2) ::
     >  T, nes

C----------------------------------------------------------------------|
C     Set plasma parameters
C----------------------------------------------------------------------|
      ne = 3e19_dp
      Te = 4e3_dp
      p  = 1._dp
      E  = 1._dp
      T(:)  = (/ 1e3_dp, 2e3_dp /)
      nes(:) = (/ 1e19_dp, 2e19_dp /)

C----------------------------------------------------------------------|
C     Test if calculation of plasma quantities working
C----------------------------------------------------------------------|
        ! ----- Coulomb logarithms ------------------------------------|
      write(*,'(a,f8.4)') 'ln Lambda_0:  ', ln_Lambda_0(ne, Te)
      write(*,'(a,10f8.4)') 'ln Lambda_0:  ', ln_Lambda_0(nes, T)
      write(*,'(a,f8.4)') 'ln Lambda_c:  ', ln_Lambda_c(ne, Te)
      write(*,'(a,f8.4)') 'ln Lambda_ee: ', ln_Lambda_ee(ne, Te, p)
      write(*,'(a,f8.4)') 'ln Lambda_ei: ', ln_Lambda_ei(ne, Te, p)

        ! ----- Thermal velocity and collision frequency --------------|
      write(*,'(a,e12.4)') 'v_th:  ', v_th(Te)
      write(*,'(a,e12.4)') 'nu_ee: ', nu_ee(ne, Te)
      write(*,'(a,e12.4)') 'v_c:   ', v_c(ne, Te, E)

        ! ----- Characteristic electric fields ------------------------|
      write(*,'(a,f8.4)') 'E_c: ', E_c(ne, Te)
      write(*,'(a,f8.4)') 'E_D: ', E_D(ne, Te)

      tmp = hot_tail_density(3.e-4_dp, 1.e-4_dp, 1._dp, 
     >  10._dp, 8.e3_dp, ne_i=3.e19_dp, ne_f=14.e19_dp, ne=5.e19_dp)

      contains


      end program runawayElectronGeneration_demo
C======================================================================|
