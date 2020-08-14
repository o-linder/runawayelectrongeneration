C######################################################################|
      module calculate_hot_tail_population
C----------------------------------------------------------------------|
C     This module is used to calculate the hot-tail population.
C
C     Functions:
C       - hot_tail_density
C       - set_ne_for_hot_tail
C       - EVDF_modified_Maxwellian
C       - tau
C       - v_c
C       - v_th
C
C     Function arguments:
C       - time: time after onset of TQ (s)
C       - ne:   electron density at time `time` (m**-3)
C       - ne_i: electron density at onset of TQ (m**-3)
C       - ne_f: electron density at end of TQ (m**-3)
C       - T:    electron temperature at time `time` (eV)
C       - T_i:  electron temperature at onset of TQ (eV)
C       - T_f:  electron temperature at end of TQ (eV)
C       - Epar: parallel electric field at time `time` (V/m)
C       - v:    electron velocity (m/s)
C
C----------------------------------------------------------------------|
      use double
      use physical_constants, only :
     >  e, epsilon_0, m_e, pi
      use Coulomb_logarithms, only :
     >  ln_Lambda_0
      use collision_frequencies, only :
     >  nu_ee

      implicit none

      public ::
     >  EVDF_modified_Maxwellian,
     >  tau,
     >  v_c,
     >  v_th

      private ::
     >  set_ne_for_hot_tail,
     >
     >  dp, e, EVDF_min_val, iter_acc, m_e, pi

        ! ----- Parameters --------------------------------------------|
      integer, parameter ::
     >  n_eval_start= 10000,            ! start # of integrand evaluatns
     >  n_eval_max  = 100000000         ! max # of integrand evaluations

      real(kind=dp) ::
     >  EVDF_min_val= 1.e-100_dp,       ! to calc. up. integration bound
     >  iter_acc    = 1.0e-8_dp         ! target accuracy of iterations

      contains

C======================================================================|
      real(kind=dp) function hot_tail_density(time, t_dec, Epar, T, T_i,
     >  ne, ne_i, ne_f) result(n_ht)
C----------------------------------------------------------------------|
C     Evaluates the density of hot-tail runaway electron.
C   
C     Input:
C       time:       time after onset of TQ (s)
C       t_dec:      time scale of exponential Te decay (s)
C       Epar:       parallel electric field at time `time` (V/m)
C       T:          electron temperature at time `time` (eV)
C       T_i:        electron temperature at onset of TQ (eV)
C       **ne:       electron density at time `time` (m**-3)
C       **ne_i:     electron density at onset of TQ (m**-3)
C       **ne_f:     electron density at end of TQ (m**-3)
C     Note, that at least one value of ne, ne_i or ne_f is required!
C     ** = optional
C
C     Output:
C       n_ht:      density of hot-tail runaway electrons (m**-3)
C
C----------------------------------------------------------------------|
        ! ----- Function arguments ------------------------------------|
      real(kind=dp), intent(in) ::
     >  time, t_dec, Epar, T, T_i
      real(kind=dp), intent(in), optional ::
     >  ne, ne_i, ne_f

        ! ----- Local variables ---------------------------------------|
      integer ::
     >  i, n_eval

      real(kind=dp) ::
     >  dv, int_bound_upp, int_bound_low, n_e, n_ef, n_ei, n_ht_old, 
     >  tau_t, v, v_th0, weight

C----------------------------------------------------------------------|
C     Initialize hot-tail density
C----------------------------------------------------------------------|
      n_ht = 0._dp

C----------------------------------------------------------------------|
C     Check input for electron density
C----------------------------------------------------------------------|
      n_e  = 0._dp;     if (present(ne  )) n_e  = ne
      n_ei = 0._dp;     if (present(ne_i)) n_ei = ne_i
      n_ef = 0._dp;     if (present(ne_f)) n_ef = ne_f

      call set_ne_for_hot_tail(time, t_dec, n_e, n_ei, n_ef)

C----------------------------------------------------------------------|
C     Determine integration bounds, normalized to thermal velocity
C----------------------------------------------------------------------|
      tau_t = tau(time, t_dec, n_ei, n_ef, T_i)
      v_th0 = v_th(T_i)

      int_bound_low = v_c(n_e, T, Epar)/v_th0
      int_bound_upp =
     >  (sqrt((-log(EVDF_min_val))**3) - 3*tau_t)**(1._dp/3._dp)

      if (int_bound_upp .le. int_bound_low) return

C----------------------------------------------------------------------|
C     Perform numerical integration using Kepler's rule
C----------------------------------------------------------------------|
        ! Set densities such that loop is executed
      n_ht     = 1._dp
      n_ht_old = 2._dp

        ! Initial number of integrand evaluations
        ! Has to be even to have uneven number of points
      n_eval = n_eval_start

        ! ----- Iteratively integrate ---------------------------------|
      do while (abs((n_ht - n_ht_old)/n_ht_old) .gt. iter_acc
     >      .and. n_eval .le. n_eval_max)
        
            ! Store previous result
        n_ht_old = n_ht

            ! Integration step size
        dv = 1._dp/n_eval*(int_bound_upp - int_bound_low)

            ! Integrand is zero for v = int_bound_low
        n_ht = 0._dp
        do i=1,n_eval
                ! Weights for integrand
            if (i.lt.n_eval) then
                weight = 2._dp**(1+mod(i,2))
            else
                weight = 1._dp
            endif

            v = i*dv + int_bound_low
            n_ht = n_ht + weight*(v**2 - int_bound_low**2)
     >          * exp(-(v**3 + 3*tau_t)**(2._dp/3._dp))
        enddo
            ! Apply pre-factors
        n_ht = 4*n_ei/sqrt(pi)*dv/3*n_ht

        write(*,'(a,i2,a,e16.8)') 'Density after step ', 
     >      1+int(log(1._dp*n_eval/n_eval_start)/log(2._dp)), 
     >      ': ', n_ht

            ! Increase number of integrand evaluations
        n_eval = 2*n_eval
      enddo

      return
      end function hot_tail_density
C======================================================================|


C======================================================================|
      subroutine set_ne_for_hot_tail(t, t_dec, y, y_i, y_f)
C----------------------------------------------------------------------|
C     Ensures the electron density has the correct format to calculate 
C     the hot-tail density.
C----------------------------------------------------------------------|
        ! ----- Parameters --------------------------------------------|
      real(kind=dp), parameter ::
     >  y_def = 1.e19_dp
        
        ! ----- Subroutine arguments ----------------------------------|
      real(kind=dp), intent(in) ::
     >  t, t_dec

      real(kind=dp), intent(inout) ::
     >  y, y_i, y_f

C----------------------------------------------------------------------|
C     Check different cases      
C----------------------------------------------------------------------|
        ! Everything provided: return
      if (y.gt.0._dp .and. y_i.gt.0._dp .and. y_f.gt.0._dp) then
        return

        ! Initial and final value provided: calculate value at `t`
      else if (y.le.0._dp .and. y_i.gt.0._dp .and. y_f.gt.0._dp) then
        y = (y_i - y_f)*exp(-t/t_dec) + y_f

        ! Initial and value at `t` provided: calculate final value
      else if (y.gt.0._dp .and. y_i.gt.0._dp .and. y_f.le.0._dp) then
        y_f = (y - y_i*exp(-t/t_dec))/(1._dp - exp(-t/t_dec))

        ! Value at `t` and final value provided: calculate initial value
      else if (y.gt.0._dp .and. y_i.le.0._dp .and. y_f.gt.0._dp) then
        y_i = (y - y_f*(1._dp - exp(-t/t_dec)))*exp(t/t_dec)

        ! Only one value provided: Set for all values
      else if (y.gt.0._dp .or. y_i.gt.0._dp .or. y_f.gt.0._dp) then
        y   = max(y, y_i, y_f)
        y_i = y
        y_f = y

        ! Nothing provided: WARNING
      else if (y.le.0._dp .and. y_i.le.0._dp .and. y_f.le.0._dp) then
        write(*,'(a,e12.4,a)') 'WARNING: No values for electron ' //
     >      'density provided for calculation of hot-tail density! ' //
     >      'Will use default density of ne=', y_def,' m**-3.'
        y   = y_def
        y_i = y_def
        y_f = y_def

      endif

      return
      end subroutine set_ne_for_hot_tail
C======================================================================|


C======================================================================|
      real(kind=dp) function EVDF_modified_Maxwellian(v, ne, v_th, tau)
     >  result(distr_func)
C----------------------------------------------------------------------|
C     Calculates the value of the Maxwellian electron velocity 
C     distribution function at velocity `v` for electron density `ne`, 
C     thermal velocity `v_th` and parameter `tau`, where the argument
C     of the exponential has been modified. Reduces to a standard
C     Maxwellian for tau=0. For details, see H.M. Smith and E. 
C     Verwichte. Phys. Plasmas 15, 072502 (2008), especially eq. (9).
C
C     Input
C       v:          electron velocity where distribution function is 
C                   to be evaluated (m/s)
C       ne:         electron density (m**-3)
C       v_th:       thermal electron velocity (m/s)
C       tau:        parameter of temporal delay
C
C     Output:
C       distr_func: magnitude of the electron distribution function
C
C----------------------------------------------------------------------|
          ! ----- Function arguments ----------------------------------|
      real(kind=dp), intent(in) ::
     >  v, ne, v_th, tau

C----------------------------------------------------------------------|
      distr_func = ne*exp(-((v/v_th)**3 + 3*tau)**(2._dp/3._dp))
     >  /(sqrt(pi)*v_th)**3

      return
      end function EVDF_modified_Maxwellian
C======================================================================|


C======================================================================|
      real(kind=dp) function tau(time, t_char, ne_i, ne_f, T_i)
C----------------------------------------------------------------------|
C     Calculates the parameter tau of temporal delay. For details, see 
C     H.M. Smith and E. Verwichte. Phys. Plasmas 15, 072502 (2008).
C
C     Input:
C       time:       time since onset of TQ (s)
C       ne:         electron density at time `time` (m**-3)
C       t_char:     characteristic time scale of TQ (s)
C       ne_i:       electron density at onset of TQ (m**-3)
C       ne_f:       electron density at end of TQ (m**-3)
C       T_i:        electron temperature at onset of TQ (eV)
C
C     Output:
C       tau:        parameter of temporal delay
C
C----------------------------------------------------------------------|
          ! ----- Function arguments ----------------------------------|
      real(kind=dp), intent(in) ::
     >  time, t_char, ne_i, ne_f, T_i

C----------------------------------------------------------------------|
      if (time .lt. 2*t_char) then
        tau = time**2/4/t_char
      else if (time .ge. 2*t_char) then
        tau = time - t_char
      endif

      tau = tau*nu_ee(ne_i, T_i)*ne_f/ne_i

      return
      end function tau
C======================================================================|


C======================================================================|
      real(kind=dp) function v_c(ne, T, Epar)
C----------------------------------------------------------------------|
C     Calculates the critical velocity for electron runaway. See e.g.
C     H.M. Smith and Epar. Verwichte. Phys. Plasmas 15, 072502 (2008).
C   
C     Input:
C       ne:         electron density (m**-3)
C       T:          electron temperature (eV)
C       Epar:       parallel electric field (V/m)
C
C     Output:
C       v_th:       thermal electron velocity (m/s)
C
C----------------------------------------------------------------------|
          ! ----- Function arguments ----------------------------------|
      real(kind=dp), intent(in) ::
     >  ne, T, Epar

C----------------------------------------------------------------------|
      v_c = sqrt(e**3/(4*pi*epsilon_0**2*m_e)
     >  *ne*ln_Lambda_0(ne, T)/Epar)

      return
      end function v_c
C======================================================================|


C======================================================================|
      real(kind=dp) function v_th(T)
C----------------------------------------------------------------------|
C     Calculates the thermal electron velocity v_th = sqrt(2*T/m_e).
C   
C     Input:
C       T:          electron temperature (eV)
C
C     Output:
C       v_th:       thermal electron velocity (m/s)
C
C----------------------------------------------------------------------|
          ! ----- Function arguments ----------------------------------|
      real(kind=dp), intent(in) ::
     >  T
C----------------------------------------------------------------------|
      v_th = sqrt(2*T*e/m_e)

      return
      end function v_th
C======================================================================|


      end module calculate_hot_tail_population
C######################################################################|
