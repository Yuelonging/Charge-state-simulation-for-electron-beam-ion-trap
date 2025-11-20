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

### Charge exchange (CX)
The Charge exchange refers that the ion decreases the charge state because if the interaction with the neutral gases, could be described by:

$B+A^{q+} = B^{n+}+A^{q-n+}$

In the project, only the single and double charge exchange are considered, 
and the [Mueller-Salzborn formula](https://doi.org/10.1016/0375-9601(77)90672-7 "The charge exhcange cross section calculation") are used to estimate the cross section

###
## ElectronBeamClass
