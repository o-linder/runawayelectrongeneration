#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# -----------------------------------------------------------------------------|
#   Header
# -----------------------------------------------------------------------------|
from    matplotlib import rc
import  matplotlib.pyplot as plt
#from    modules.dataTools import max_along_axis
import  numpy as np
import  scipy.constants as physConst
import  scipy.integrate as integrate
from    scipy.optimize import curve_fit

plt.ion()
rc('text', usetex=True)
rc('font', size=10, family='serif')

# -----------------------------------------------------------------------------|
#   Class to calculate hot-tail population
# -----------------------------------------------------------------------------|
class hot_tail_generation:
        # minimum value for exponential during integration
    __exp_min_val = 1e-100      

        # factor for converting number to current density 
    __j_conv = physConst.e*physConst.c

        # maximum number of iteration during quadrature
    __quad_max_iter = 1000 

    # -------------------------------------------------------------------------|
    def __init__(self, t, E, t_dec=None, t_del=0, ne=None, ne_i=None, 
            ne_f=None, Te=None, Te_i=None, Te_f=None, calc_evolution=True):
            # ----- Can hot-tail calculations be performed --------------------|
        self.calc_possible = True

            # ----- Time array, delay and decay -------------------------------|
        self.t = np.atleast_1d(t)
        self.t_dec = t_dec
        self.t_del = t_del

            # ----- Electric field --------------------------------------------|
        self.E = np.abs(np.atleast_1d(E))
        if self.E.size == 1:
            self.E *= np.ones(self.t.shape)

            # ----- Electron temperature --------------------------------------|
        if Te is not None:
            self.Te = np.atleast_1d(Te)

            if self.Te.size == 1:
                self.Te_i = Te_i
                self.Te_f = Te_f

            elif self.Te.size > 1:
                self.Te_i = self.Te[0]
                self.Te_f = self.Te[-1]
                if self.t_dec is None:
                    print('Decay time not provided. Trying to perform a fit.')
                    self.t_dec =self.fit_exponential(self.t, self.Te)[0] 

        elif np.all(np.array([Te_i, Te_f, self.t_dec]) != None):
            self.Te_i = Te_i
            self.Te_f = Te_f
            self.Te = self.Te_f + (self.Te_i - self.Te_f)\
                    *np.exp(-self.t/self.t_dec)
        else:
            self.calc_possible = False
            print('Cannot set electron temperature.')

            # ----- Electron density ------------------------------------------|
        self.set_electron_density(ne=ne, ne_i=ne_i, ne_f=ne_f)

            # ----- Additional quantities -------------------------------------|
        self.nu_0 = np.zeros(self.t.shape)
        self.v_T0 = np.zeros(self.t.shape)
        self.v_c = np.zeros(self.t.shape)
        self.tau = np.zeros(self.t.shape)
        self.calc_additional_quantities()

            # ----- Calculate evolution of the hot-tail population ------------|
        self.n_hot = np.zeros(self.t.shape)
        self.j_hot = np.zeros(self.t.shape)
        if calc_evolution:
            self.calc_evolution()

        # ----- end method __init__ -------------------------------------------|

    # -------------------------------------------------------------------------|
    def calc_evolution(self, assume_single_max=False, increasing_only=True):
        """Calculates the evolution of the hot-tail population. If the switch 
        `assume_single_max` is set, the calculation is stopped as soon as the 
        first maximum is encountered.
        """
        self.n_hot = np.zeros(self.t.shape)

            # Check if hot-tail calculation possible
        if not self.calc_possible:
            print('Calculation of hot-tail population not possible. Abort.')
            return

            # ----- Evolve hot-tail population --------------------------------|
        for i in range(self.t.size):
            if self.t[i] < self.t_del: continue

                # ----- Determine integration limits --------------------------|
                # Between v_c and where exponential drops below a value of 
                # `__exp_min_val`
            int_lim = ( self.v_c[i], 
                        ((-np.log(self.__exp_min_val))**(3/2)-3*self.tau[i])**(1/3)\
                                *self.v_T0)
            if int_lim[1]/self.v_c[i] < 1 or np.isnan(int_lim[1]): continue

                # ----- Hot-tail population at `t[i]` -------------------------|
            self.n_hot[i] = 4*self.ne_i/(np.sqrt(np.pi)*self.v_T0**3) \
                    *integrate.quadrature(
                            lambda v: np.exp(-((v/self.v_T0)**3 + 3*self.tau[i])**(2/3)) \
                                    *(v**2 - self.v_c[i]**2),
                            *int_lim, maxiter=self.__quad_max_iter)[0] 
        
                        # stop calculation if maximum has been reached
            if assume_single_max and i > 0 and self.n_hot[i] < self.n_hot[i-1]:
                break

            # ----- Final hot-tail density does not decay ---------------------|
            # This assumes, that electrons with velocities exceeding the 
            # critical velocity do not equilibriate through collisions since
            # they experience net acceleration by the applied electric field.
#        if increasing_only:
#            __ = max_along_axis(self.n_hot)

            # ----- Calculate hot-tail carried current ------------------------|
            # This assumes j_hot = e c n_hot
        self.j_hot = self.__j_conv * self.n_hot

        # ----- end method calc_evolution -------------------------------------|


    # -------------------------------------------------------------------------|
    #   Setup electron temperature and density profiles
    # -------------------------------------------------------------------------|
    def set_electron_density(self, ne=None, ne_i=None, ne_f=None):
        """Function to set the electron density evolution.
        """
        if ne is not None:
            self.ne = np.atleast_1d(ne)

            if self.ne.size == 1:
                self.ne_i = ne_i
                self.ne_f = ne_f
            elif self.ne.size > 1:
                self.ne_i = self.ne[0]
                self.ne_f = self.ne[-1]

        elif np.all(np.array([ne_i, ne_f, self.t_dec]) != None):
            self.ne_i = ne_i
            self.ne_f = ne_f
            self.ne = self.ne_f + (self.ne_i - self.ne_f)\
                    *np.exp(-self.t/self.t_dec)

        elif ne_i is not None:
            self.ne_i = ne_i
            self.ne_f = ne_i
            self.ne = ne_i*np.ones(self.t.shape)

        else:
            self.calc_possible = False
            print('Cannot set electron density. Abort.')

        # ----- end method set_electron_density -------------------------------|


    def fit_exponential(self, x, y):
        """Fit an exponential to the data (`x`, `y`) by taking the logarim of 
        `y` and fitting a linear function to it, thus retrieve the decay time.
        """
        popt, pcov = curve_fit(self.lin_func, x, np.log(y), p0=(1e-4, 1e0))
        return popt[0], np.sqrt(pcov[0,0])
        # ----- end method fit_exponential ------------------------------------|

    # -------------------------------------------------------------------------|
    def lin_func(self, x, a, b):
        """Linear function for interpolation, yielding the negative, inverse 
        slope `a` and the offset `b`. This can be used to determine a decay 
        time for an exponentially decreasing function.
        """
        return -x/a+b
        # ----- end method lin_func -------------------------------------------|

    # -------------------------------------------------------------------------|
    #   Additional quantities necessary to determine hot-tail population   
    # -------------------------------------------------------------------------|
    def calc_additional_quantities(self):
        """Calculates additional quantities needed to evaluate the evolution 
        of the hot-tail population.
        """
        if not self.calc_possible: return

                # initial collision frequency
        self.nu_0 = self.__nu__(self.ne_i, self.Te_i)

                # initial thermal velocity
        self.v_T0 = self.__v_T__(self.Te_i)

                # critical velocity
        self.v_c = self.__v_c__(self.ne, self.Te, self.E)

                # tau
        self.tau = self.__tau__(self.t, self.t_dec, self.nu_0, 
                ne_i=self.ne_i, ne_f=self.ne_f, method='ppg')

        # ----- end method calc_additional_quantities -------------------------|

        # ---------------------------------------------------------------------|
    def __EVDF__(self, v, n, v_T, tau=0):
        """Calculates the value of the Maxwellian electron velocity 
        distribution function at velocity `v` in units of m/s for electron 
        density `n` in units of m**-3, thermal velocity `v_T` in units of m/s
        and `tau`.

        From H.M. Smith and E. Verwichte. Phys. Plasmas 15, 072502 (2008), 
        eq. (9).
        """
        return n/(np.sqrt(np.pi)*v_T)**3*np.exp(-((v/v_T)**3 + 3*tau)**(2/3))
        # ----- end method __EVDF__ -------------------------------------------|

        # ---------------------------------------------------------------------|
    def __lnLambda__(self, n, T):
        """
        Calculates Coulomb logarithm for electron-electron collisions of thermal particles of 
        density `n` in units of m**-3 and temperature `T` in units of eV.

        From J. Wesson. Tokamaks. Oxford University Press 2004, p. 727.
        """
        return 14.9 - .5*np.log(n*1e-20) + np.log(1e-3*T)
        # ----- end method __lnLambda__ ---------------------------------------|

        # ---------------------------------------------------------------------|
    def __nu__(self, n, T):
        """
        Calculates the electron-electron collision frequency for thermal particles of density
        `n` in units of m**-3 and temperature `T` in units of eV.

        From P. Helander et al., Plasma Phys. Control. Fusion 44, B247 (2002).
        """
        return n*self.__lnLambda__(n, T)/self.__v_T__(T)**3 \
            *physConst.e**4/(4*np.pi*physConst.epsilon_0**2*physConst.m_e**2)
        # ---- end method __nu__ ----------------------------------------------|

        # ---------------------------------------------------------------------|
    def __tau__(self, t, t_char, nu_0, ne_i=1, ne_f=0, method='ppg'):
        """
        Calcualtes the parameter tau for hot-tail generation using either the `method` 'ppg' from
        Geri's implementation or 'Smith' from H.M. Smith and E. Verwichte. Phys. Plasmas 15, 072502 
        (2008), eq. (17). In case of 'ppg', the characteristic time `t_char` is the exponential 
        decay time, in case of 'Smith', `t_char` is the time delay.
        """
            # ----- Check input -----------------------------------------------|
            # Eliminates the need of providing initial and final electron 
            # density if this quantity does not change throughout the 
            #  temperature decay.
        if ne_f == 0:
            ne_f = ne_i

            # ----- Calculate quantity tau ------------------------------------|
        tau = np.empty(t.shape)
        if method=='ppg':
            tau[t <  2*t_char] = t[t <  2*t_char]**2/4/t_char
            tau[t >= 2*t_char] = t[t >= 2*t_char] - t_char

        elif method=='Smith':
            tau[t <= t_char] = 0.
            tau[t >  t_char] = t[t > t_char] - t_char

        return tau*nu_0*ne_f/ne_i
        # ----- end method __tau__ --------------------------------------------|

        # ---------------------------------------------------------------------|
    def __v_c__(self, n, T, E):
        """
        Calculates critical velocity for electron runaway with electron density `n` in units of 
        m**-3, electron temperature `T` in units of eV and external electric field `E` in units
        of V/m.

        From H.M. Smith and E. Verwichte. Phys. Plasmas 15, 072502 (2008).
        """
        return np.sqrt(n*physConst.e**3*self.__lnLambda__(n, T)) \
                /np.sqrt((4*np.pi*physConst.epsilon_0**2*physConst.m_e*E))

        # ---------------------------------------------------------------------|
    def __v_T__(self, T):
        """
        Calculates electron thermal velocity at temperature `T`, with `T` in units of eV.
        """
        return np.sqrt(2*T*physConst.e/physConst.m_e)
        # ----- end method __v_T__ --------------------------------------------|

    # -------------------------------------------------------------------------|
    #   Plot the evolution of key quantities, being the 
    # -------------------------------------------------------------------------|
    def plot_evolution(self):
        """
        Plot the evolution of the hot-tail population and associated quantities.
        """
        fig, ax = plt.subplots(3, 2, figsize=(7,6))
        ax = fig.axes
            
        ax[0].plot(self.t, 1e-16*self.n_hot, c='k')
        ax[0].set_title(r'Hot-tail population')
        ax[0].set_ylabel(r'$n_{\rm hot}$~(10$^{16}$~ m$^{-3}$)')
        ax_t = ax[0].twinx()
        ax_t.plot(self.t, 1e-6*self.j_hot, c='k')
        ax_t.set_ylabel(r'$j_{\rm hot}$~(MA/m$^2$)')
        ax_t.set_ylim(bottom=0)

        ax[1].plot(self.t, self.Te, c='k')
        ax[1].semilogy()
        ax[1].set_title('Electron temperature')
        ax[1].set_ylabel(r'$T_{\rm e}$~(eV)')
        ax[1].set_ylim(bottom=1)

        ax[2].plot(self.t, self.v_c/self.v_T0, c='k')
        ax[2].set_title('Critical velocity')
        ax[2].set_ylabel(r'$v_{\rm c}/v_{T_0}$')

        ax[3].plot(self.t, 1e-19*self.ne, c='k')
        ax[3].set_title('Electron density')
        ax[3].set_ylabel(r'$n_{\rm e}$~(10$^{19}$~m$^{-3}$)')

        ax[4].plot(self.t, self.tau, c='k')
        ax[4].set_title(r'$\tau$')
        ax[4].set_ylabel(r'$\tau$')

        ax[5].plot(self.t, self.E, c='k')
        ax[5].set_title('Electric field')
        ax[5].set_ylabel(r'$E$~(V/m)')
        
        for i, a in enumerate(ax):
            a.set_xlabel(r'$t~({\rm s})$')
            a.set_xlim((self.t[0], self.t[-1]))

            if i != 1:
                a.set_ylim(bottom=0)

        plt.tight_layout()

        return fig
        # ----- end method plot_evolution -------------------------------------|

# -----------------------------------------------------------------------------|
#   Function to demonstrate hot-tail population evolution
# -----------------------------------------------------------------------------|
def demo():
    t = np.arange(0, 2.e-3 + 5.e-6, 5.e-6)
    E = 1. + (0.01 - 1.)*np.exp(-t/5.e-4)
    ht = hot_tail_generation(t, E, t_del=0, t_dec=1.5e-4, 
            ne_i=3.e19, ne_f=15.e19, Te_i=7.e3, Te_f=10, calc_evolution=False)
    ht.calc_evolution(assume_single_max=False, increasing_only=False)

#    ht.plot_evolution()

    return ht
    # ----- end function demo -------------------------------------------------|

# -----------------------------------------------------------------------------|
#   Run demo
# -----------------------------------------------------------------------------|
ht = demo()

np.savetxt('dat/hot_tails_python.dat', 
        np.array([ht.t, ht.n_hot, ht.ne, ht.Te, ht.E, ht.v_c/ht.v_T0, ht.tau]).T,
        fmt='%19.12e',
        header=   'Time (s)          ' + \
                '  n_hot (m**-3)     ' + \
                '  n_e (m**-3)       ' + \
                '  T_e (ev)          ' + \
                '  E_par (V/m)       ' + \
                '  v_c (v_th0)       ' + \
                '  tau',
        )

# ----- end script hot_tails.py -----------------------------------------------|
