C======================================================================|
      program hottailpopulation_demo
C----------------------------------------------------------------------|
C     Demonstration of calculating the hot-tail population.
C----------------------------------------------------------------------|
      use double
      use runawayelectrongeneration

      implicit none

        ! ----- Parameters --------------------------------------------|
      integer, parameter ::
     >  n_time = 401

      real(kind=dp), parameter ::
     >  time_step = .5e-5_dp

        ! ----- Local variables ---------------------------------------|
      integer ::
     >  i, un

      integer, dimension(2) ::
     >  dings

      real(kind=dp) ::
     >  ne_i, ne_f, t_dec, Te_i, Te_f, tmp

      real(kind=dp), dimension(n_time) ::
     >  Epar, ne, Te, time

C----------------------------------------------------------------------|
C     Set plasma parameters
C----------------------------------------------------------------------|
      time  = (/((i-1)*time_step, i=1,n_time)/)
      t_dec = 1.5e-4_dp
      Te_i  = 7.e3_dp
      Te_f  = 10._dp
      ne_i  = 3.e19_dp
      ne_f  = 15.e19_dp

      Te    = (Te_i - Te_f)*exp(-time/t_dec) + Te_f
      ne    = (ne_i - ne_f)*exp(-time/t_dec) + ne_f
      Epar  = (0.01_dp - 1._dp)*exp(-time/5.e-4_dp) + 1._dp

C----------------------------------------------------------------------|
C     Calculate hot-tail density for demo parameters
C----------------------------------------------------------------------|
      un=10
      open(unit=un, file='hot_tails.dat', action='write')

        ! Header
      write(un,'(a)',advance='no') '# Time (s)          ' 
      write(un,'(a)',advance='no') '  n_hot (m**-3)     '
      write(un,'(a)',advance='no') '  n_e (m**-3)       '
      write(un,'(a)',advance='no') '  T_e (ev)          '
      write(un,'(a)',advance='no') '  E_par (V/m)       '
      write(un,'(a)',advance='no') '  v_c (v_th0)       '
      write(un,'(a)')              '  tau               '

      do i=1,n_time
        tmp = hot_tail_density(time(i), t_dec, Epar(i), T=Te(i), 
     >      T_i=Te_i, ne=ne(i), ne_i=ne_i, ne_f=ne_f, 
     >      store_result=.true., un=un) 
      enddo
      close(un)

      contains

      end program hottailpopulation_demo
C======================================================================|
