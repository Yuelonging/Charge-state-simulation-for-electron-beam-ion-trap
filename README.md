#  Project struce
This program consists of three classes: ElectronBeam, Element, Plasma. By pre-defining the parameters, 
you could simulate the charge state distribution in the with a given electron energy.

## The ion charge state 
In the following chaper, I will introduce the physical progress determing the charge state that hase been considered.

### Electron impact ionization (EII) 
The EII progress could be described by the following progrees:

$e^-(E)+A^{q+}=e^-(E_1)+e^-(E_2)+A^{q+1+}$

In this project, the [lotz experience formula](https://doi.org/10.1007/BF01325928 "The electron impact ionization cross section calculation") 
was used to calculate the EII cross section $(cm^{-2})$

### Radiative recombination (RR)
The RR progress could be described by:

$e^-+A^{q+}=A^{q-1+}+\gamma$

[The Radiative combination cross section](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.27.2913 "The radiative combination cross section calculation")
is calculated by the PhysRevA.27.2913

### Ion Escape
The ion escape refers that the ion escaping from the confinement provided by the electron space charge effect and the potential well applied on the drift tubes.
The escape rate just takes the solution of the [planck-fock equation with the square potential well](). The escape rate is related to the ion temperature,
and both radial and axial escaping are considered.

### Charge exchange (CX)
The Charge exchange refers that the ion decreases the charge state because if the interaction with the neutral gases, could be described by:

$B+A^{q+} = B^{n+}+A^{q-n+}$

In the project, only the single and double charge exchange are considered, 
and the [Mueller-Salzborn formula](https://doi.org/10.1016/0375-9601(77)90672-7 "The charge exhcange cross section calculation") are used to estimate the cross section

## The ion temperature
In the following chapter, I will introduce the physcial progress that determines the ion temperature.

### The electron heating
The electron heating (Spitzer heating) comes from the collisional energy exchange between the electrons and ions:

$\frac{dE_i}{dt}=\frac{4}{3}\nu_{ei}N_i\frac{m_e}{M_i}E_e$

## ElectronBeamClass
