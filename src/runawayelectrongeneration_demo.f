C======================================================================|
      program runawayElectronGeneration_demo

      use double
      use runawayelectrongeneration

      implicit none

      real(kind=dp) ::
     >  ne, p, Te

C----------------------------------------------------------------------|
C     Set plasma parameters
C----------------------------------------------------------------------|
      ne = 3e19_dp
      Te = 4e3_dp
      p  = 1._dp

C----------------------------------------------------------------------|
C     Test if calculation of plasma quantities working
C----------------------------------------------------------------------|
        ! ----- Coulomb logarithms ------------------------------------|
      write(*,'(a,f8.4)') 'ln Lambda_0:  ', ln_Lambda_0(ne, Te)
      write(*,'(a,f8.4)') 'ln Lambda_c:  ', ln_Lambda_c(ne, Te)
      write(*,'(a,f8.4)') 'ln Lambda_ee: ', ln_Lambda_ee(ne, Te, p)
      write(*,'(a,f8.4)') 'ln Lambda_ei: ', ln_Lambda_ei(ne, Te, p)

        ! ----- Thermal collision frequency ---------------------------|
      write(*,'(a,e12.4)') 'nu_ee: ', nu_ee(ne, Te)

        ! ----- Characteristic electric fields ------------------------|
      write(*,'(a,f8.4)') 'E_c: ', E_c(ne, Te)
      write(*,'(a,f8.4)') 'E_D: ', E_D(ne, Te)

      end program runawayElectronGeneration_demo
C======================================================================|
