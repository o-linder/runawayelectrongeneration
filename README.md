# Runaway Electron Generation 
This project provides Fortran procedures for the calculation of runaway electron growth rates due to momentum-space diffusion (Dreicer generation) and knock-on collisions (avalanche generation), as well as for the calculation of a hot-tail runaway electron population.

## Getting Started
Compile the .f-files located in `src` in the following order to generate the project object and module files:

1. `double.f`
2. `physical_contants.f`
3. `file_io.f`
4. `Coulomb_logarithms.f`
5. `collision_frequencies.f`
6. `electric_fields.f`
7. `calculate_hot_tail_population.f`
8. `calculate_Dreicer_growthrate.f`
9. `calculate_avalanche_growthrate.f`
10. `runawayelectrongeneartion.f`

An example is given in the makefile (`make project`/`make all`). 

A demonstration program for the calculation of the hot tail population is provided in `src/hot_tail_demo.f` and may be compiled using the makefile (`make hot_tail_demo`/`make all`).

*Note, project developed using the Intel Fortran compiler `ifort 14.0`. Building project with other compilers not tested.*

## Using the Runaway Electron Generation Module
The module `runawayelectrongeneration` can be used for the calculation of the

1.  Hot-tail runaway electron population,
2.  Dreicer runaway electron growth rate,
3.  Avalanche runaway electron growth rate.

Throughout the project, SI units are employed, except in the case of temperatures, which are used units of eV.

### Function Arguments
The following arguments are used by the functions in this module:

* `real(kind=dp) B`: magnetic field (T)
* `real(kind=dp) Epar`: parallel electric field (V/m)
* `real(kind=dp) eps`: inverse tokamak aspect ratio 
* `real(kind=dp) ne`: electron density at time `T` (m\*\*-3)
* `real(kind=dp) ne_f`: electron density at end of thermal quench (m\*\*-3)
* `real(kind=dp) ne_i`: electron density at onset of thermal quench (m\*\*-3)
* `real(kind=dp) nj(nspecies)`: density of each ion species (m\*\*-3)
* `integer nspecies`: number of ion species
* `real(kind=dp) T`: electron temperature (at time `T`) (eV)
* `real(kind=dp) t_dec`: time scale of electron temperature decay (s)
* `real(kind=dp) T_f`: electron temperature at end of thermal quench (eV)
* `real(kind=dp) T_i`: electron temperature at onset of thermal quench (eV)
* `real(kind=dp) time': time after onset of thermal quench (s)
* `real(kind=dp) v`: electron velocity (m/s)
* `integer Z(nspecies)`: atomic number of each ion species
* `integer Z0(nspecies)`: net charge of each ion species

### Hot-tail Runaway Electron Population
The function 
```fortran
      real(kind=dp) function hot_tail_density(time, t_dec, Epar, T, T_i,
     >  T_f, ne, ne_i, ne_f, verbose, store_result, un) result(n_ht)
```
calculates the hot-tail density `n_ht` in units of m\*\*-3 according to Eq. (18) in [H.M. Smith and E. Verwichte. Hot tail runaway electron generation in tokamak disruptions. *Phys. Plasmas* **15**, 072502 (2008)](https://doi.org/10.1063/1.2949692). 

Notes on the function arguments:
* The time `time` is provided with respect to the onset of the thermal quench. 
* For the electron temperatures `T`, `T_i` and `T_f`, only `T` and `T_i` are required to calculate the hot-tail population. However, when providing `T_f` and only one of the temperatures `T` and`T_i`, the remaining quantity is estimated from the equation `T(time) = (T_i - T_f) * exp(-time/t_dec) + T_i`.
* For the electron densities `ne`, `ne_i` and `ne_f`, all three quantities should ideally be provided. However, calculations can also be performed with only two quantities provided by estimating the remaining quantity from `ne(time) = (ne_i - ne_f) * exp(-time/t_dec) + ne_i`. In the case where only one density is provided, a density constant in time is assumed (not recommended).
* Optional: The flag `logical verbose` allows to print calculation parameters and results to screen.
* Optional: The flag `logical store_result` allows to store calculation parameters and results in an opened, writeable file with unit number `integer un`.
* Optional: The unit number `integer un` of an opened, writeable file has to be provided if results are to be stored. 

### Dreicer Runaway Electron Growth Rate: Classic Formula
The function 
```fortran
      real(kind=dp) function Dreicer_growthrate_classic(Z, Z0, nj, 
     >  nSpecies, T, Epar) result (Gamma_D) 
```
calculates the runaway electron growth rate `Gamma_D` in units of 1/s due to Dreicer generation according to Eqs. (62)-(64) in [J.W. Connor and R.J. Hastie. Relativistic limitations on runaway electrons, *Nucl. Fusion* **15** (1975), 415](https://doi.org/10.1088/0029-5515/15/3/007). Note, that the full expression is evaluated. Neither the limit of strong electric fields, equation (66), or the non-relativistic limit, equation (67), is included, as these are deemed irrelevant for application.

### Dreicer Runaway Electron Growth Rate: Neural Network
The function
```fortran
      real(kind=dp) function Dreicer_growthrate_CODE_neural_network(
     >   Z, Z0, nj, nSpecies, T, Epar) result (Gamma_D)
```
calculates the runaway electron growth rate `Gamma_D` in units of 1/s due to Dreicer generation employing a neural network presented in [L. Hesslow, L. Unnerfelt, O. Vallhagen, O. Embreus, M. Hoppe, G. Papp and T. Fulop. Evaluation of the Dreicer runaway growth rate in the presence of high-Z impurities using a neural network, *J. Plasma Phys* **85** (2019), 475850601](https://doi.org/10.1017/S0022377819000874). The original MATLAB implementation of the neural network can be found at [https://github.com/unnerfelt/dreicer-nn](https://github.com/unnerfelt/dreicer-nn).

### Avalanche Runaway Electron Generation Rate: Classic Formula
The function
```fortran
      real(kind=dp) function avalanche_growthrate_classic(Z, Z0, nj, 
     >  nSpecies, T, Epar, eps) result (Gamma_av)
```
calculates the runaway electron growth rate `Gamma_av` in units of 1/s due to avalanche generation according to Eq. (18) in [M.N. Rosenbluth and S.V. Putvinski. Theory for avalanche of runaway electrons in tokamaks, *Nucl. Fusion* **37** (1997), 1355](https://doi.org/10.1088/0029-5515/37/10/I03).

### Avalanche Runaway Electron Generation Rate: Including Partially Ionized Impurities
The function
```fortran
      real(kind=dp) function avalanche_growthrate(Z, Z0, nj, nSpecies, 
     >  T, Epar, B) result (Gamma_av) 
```
calculates the runaway electron growth rate `Gamma_av` in units of 1/s due to avalanche generation considering the impact of partially ionized impurities according to Eq. (14) in [L. Hesslow, O. Embreus, O. Vallhagen and T. Fulop. Influence of massive material injection on avalanche runaway generation during tokamak disruptions, *Nucl. Fusion* **59** (2019), 084004](https://doi.org/10.1088/1741-4326/ab26c2).

### Further functions and Subroutines
The module contains additional public functions and subroutines for the calculation of quantities required in aforementioned functions for the calculation of the hot-tail runaway electron population and the runaway electron growth rates. These are listed in the following sections. If not in any of the references mentioned above, these additional quantities are defined in the following publications.

1. [L. Hesslow, O. Embreus, M. Hoppe, T.C. DuBois, G. Papp, M. Rahm and T. Fulop. Generalized collision operator for fast electrons interacting with partially ionized impurities, *J. Plasma Phys.* **84** (2018), 905840605](https://doi.org/10.1017/S0022377818001113)
2. [L. Hesslow, O. Embreus, G.J. Wilkie, G. Papp and T. Fulop. Effect of partially ionized impurities and radiation on the effective critical electric field for runaway generation, *Plasma Phys. Control. Fusion* **60** (2018), 074010](https://doi.org/10.1088/1361-6587/aac33e)
3. [P. Helander, L.-G. Eriksson and F. Andersson. Runaway acceleration during magnetic reconnection in tokamaks, *Plasma Phys. Control. Fusion* **44** (2002), B247](https://doi.org/10.1088/0741-3335/44/12B/318)

#### `function E_ceff_over_Ec`
```fortran
C======================================================================|
      real(kind=dp) function E_ceff_over_E_c(Z, Z0, nj, nSpecies,
     >  T, B)
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
C
C----------------------------------------------------------------------|
```

#### `function E_c`
```fortran
C======================================================================|
      real(kind=dp) function E_c(ne, T)
C----------------------------------------------------------------------|
C     Calculates the critical electric field E_c:
C       E_c = me*c/e/tauc = ne*e**3*ln Lambda_c/(4*pi*eps0**2*me*c**2).
C
C     Input:
C       ne  : free electron density (m**-3)
C       T   : electron temperature (eV)
C
C     Output:
C       E_c : critical electric field (V)
C
C----------------------------------------------------------------------|
```

#### `function E_D`
```fortran
C======================================================================|
      real(kind=dp) function E_D(ne, T)
C----------------------------------------------------------------------|
C     Calculates the Dreicer electric field E_D:
C       E_D = ne*e**3*ln Lambda/(4*pi*eps0**2*T).
C
C     Input:
C       ne  : free electron density (m**-3)
C       T   : electron temperature (eV)
C
C     Output:
C       E_D : Dreicer electric field (V)
C
C----------------------------------------------------------------------|
```

#### `function tau`
```fortran
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
```

#### `function EVDF_modified_Maxwellian`
```fortran
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
```

#### `function p_star`
```fortran
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
```

#### `function nuBar_D`
```fortran
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
```

#### `subroutine calc_nuBar_D01`
```fortran
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
```

#### `function nuBar_S`
```fortran
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
```

#### `subroutine calc_nuBar_S01`
```fortran
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
```

#### `function nu_ee`
```fortran
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
```

#### `function v_c`
```fortran
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
```

#### `function v_th`
```fortran
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
```

#### `function ln_Lambda_0`
```fortran
C======================================================================|
      real(kind=dp) function ln_Lambda_0(ne, T)
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
C
C----------------------------------------------------------------------|
```

#### `function ln_Lambda_c`
```fortran
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
C
C----------------------------------------------------------------------|
```

#### `function ln_Lambda_ee`
```fortran
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
C
C----------------------------------------------------------------------|
```

#### `function ln_Lambda_ei`
```fortran
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
C
C----------------------------------------------------------------------|
```

