# Runaway Electron Generation 
This projects implements functions for the calculation of runaway electron growth rates from analytical expressions in Fortran.

## Getting Started
Simply compile using your favourite Fortran compiler

1. the module contained in `double.f`, 
2. the module contained in `runawayelectron.f` (which requires the module `double`).

The `runawayelectrongeneration` module can now be used in a Fortran project of yours.

## Using the Runaway Electron Generation Module
The runaway electron growth rate may be calculated for 

1.  Dreicer generation,
2.  avalanche generation.

For each generation process, two public functions are available, either omitting or including the impact of partially ionized impurities on runaway electron generation. All functions take the following input quantities

* `integer Z(ns)`: atomic number of each ion species
* `integer Z0(ns)`: net charge of each ion species
* `double n(ns)`: density of each ion species (m\*\*-3)
* `integer ns`: number of ion species
* `double T`: temperature (eV)
* `double Epar`: parallel electric field (V/m)

For the calculation of avalanche growth rates, additional quantities are necessary, depending on the method used. All functions return the runaway electron growth rate in units of 1/s.

### Dreicer Generation
Primary runaway electron generation can be calculated using one of two functions.

#### Analytical formulae: `calc_Gamma_D(Z,Z0,n,ns,T,Epar)`
This function calculates the growth rate according to equations (62-64) in [J.W. Connor and R.J. Hastie. Relativistic limitations on runaway electrons, *Nucl. Fusion* **15** (1975), 415](https://doi.org/10.1088/0029-5515/15/3/007). Note, that the full expression is evaluated. Neither the limit of strong electric fields, equation (66), or the non-relativistic limit, equation (67), is included, as these are deemed irrelevant for application.

#### Neural network: `calc_Gamma_D_nn(Z,Z0,n,ns,T,Epar)`
This function calculates the growth rate by evaluating a neural network, described in [L. Hesslow, L. Unnerfelt, O. Vallhagen, O. Embreus, M. Hoppe, G. Papp and T. Fulop. Evaluation of the Dreicer runaway growth rate in the presence of high-Z impurities using a neural network, Submitted for publication in *J. Plasma Phys*](https://arxiv.org/pdf/1910.00356). The original MATLAB implementation of the neural network can be found at [https://github.com/unnerfelt/dreicer-nn](https://github.com/unnerfelt/dreicer-nn).

### Avalanche Generation
Secondary runaway generation due to the avalanche mechanism can be calculated using one of two function.

#### Analytical formulae: `calc_Gamma_av_RP(Z,Z0,n,ns,T,Epar,eps)`
This function calculates the growth rate from equation (18) in [M.N. Rosenbluth and S.V. Putvinski. Theory for avalanche of runaway electrons in tokamaks, *Nucl. Fusion* **37** (1997), 1355](https://doi.org/10.1088/0029-5515/37/10/I03). As additional input quantity, the local aspect ratio `double eps` is required.

#### Analytical formulae with partially ionized impurities: `calc_Gamma_av(Z,Z0,n,ns,T,Epar,B)`
This function calculates the growth rate from equation (14) in [L. Hesslow, O. Embreus, O. Vallhagen and T. Fulop. Influence of massive material injection on avalanche runaway generation during tokamak disruptions, *Nucl. Fusion* **59** (2019), 084004](https://doi.org/10.1088/1741-4326/ab26c2). As additional input quantity, the local magnetic field `double B` is required.

### Further functions and subroutines
The module contains additional public functions and subroutines for the calculation of quantities required in aforementioned functions for calculation of runaway electron growth rates, being

* characteristic electric fields, being 
    - the normalized effective critical electric field `calc_EceffOverEc(Z,Z0,n,ns,T,B)` (see [https://github.com/hesslow/Eceff](https://github.com/hesslow/Eceff) for the original implementation in MATLAB), 
    - the critical electric field `calc_Ec(ne,T)`, 
    - the Dreicer electric field `calc_ED(ne,T)`,
* critical momentum for electron runaway `calc_pstar(Z,Z0,n,ns,T,Epar,B)`,
* collision frequencies, i.e. 
    - the normalized electron deflection frequency `calc_nuBarD(Z,Z0,n,ns,T,p)`, 
    - the first two terms of the expansion of the normalized electron deflection frequency `subroutine calc_nuBarD01(Z,Z0,n,ns,T,nuBarD0,nuBarD1)`, 
    - the normalized electron slowing down frequency `calc_nuBarS(Z,Z0,n,ns,T)`, 
    - the first two terms of the expansion of the normalized electron slowing down frequency `subroutine calc_nuBarS01(Z,Z0,n,ns,T,nuBarS0,nuBarS1)`, 
    - the thermal electron-electron collision frequency `calc_nuee(ne,T)`,
* Coulomb logarithms for 
    - thermal electron-electron collisions `calc_lnLambda0(ne,T)`, 
    - relativistic electron-electron collisions `calc_lnLambdac(ne,T)`, 
    - collision of non-thermal and thermal electrons `calc_lnLambdaee(ne,T,p)`,
    - collisions of non-thermal electrons and thermal ions `calc_lnLambdaei(ne,T,p)`.

The input quantity `double p` appearing in several functions describes the normalized electron momentum. Some of the quantities calculated in these functions and subroutines are defined in aforementioned references. Remaining quantites are defined in the following references:

1. [L. Hesslow, O. Embreus, M. Hoppe, T.C. DuBois, G. Papp, M. Rahm and T. Fulop. Generalized collision operator for fast electrons interacting with partially ionized impurities, *J. Plasma Phys.* **84** (2018), 905840605](https://doi.org/10.1017/S0022377818001113)
2. [L. Hesslow, O. Embreus, G.J. Wilkie, G. Papp and T. Fulop. Effect of partially ionized impurities and radiation on the effective critical electric field for runaway generation, *Plasma Phys. Control. Fusion* **60** (2018), 074010](https://doi.org/10.1088/1361-6587/aac33e)
3. [P. Helander, L.-G. Eriksson and F. Andersson. Runaway acceleration during magnetic reconnection in tokamaks, *Plasma Phys. Control. Fusion* **44** (2002), B247](https://doi.org/10.1088/0741-3335/44/12B/318)
