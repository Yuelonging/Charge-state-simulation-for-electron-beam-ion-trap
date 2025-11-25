#  Project struce
This program consists of three classes: ElectronBeam, Element, Plasma. By pre-defining the parameters, 
you could simulate the charge state distribution in the with a given electron energy.

## The ion charge state 
In the following chaper, I will introduce the physical progress determing the charge state that hase been considered.

### Electron impact ionization (EII) 
The EII progress could be described by the following progrees:

$$e^-(E)+A^{q+}=e^-(E_1)+e^-(E_2)+A^{q+1+}$$

In this project, the [lotz experience formula](https://doi.org/10.1007/BF01325928 "The electron impact ionization cross section calculation") 
was used to calculate the EII cross section $(cm^{-2})$

### Radiative recombination (RR)
The RR progress could be described by:

$$e^-+A^{q+}=A^{q-1+}+\gamma,$$

[The Radiative combination cross section](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.27.2913 "The radiative combination cross section calculation")
is calculated by the PhysRevA.27.2913

### Ion Escape
The ion escape refers that the ion escaping from the confinement provided by the electron space charge effect and the potential well applied on the drift tubes.
The escape rate just takes the solution of the [the Fokker-Plank equation with the square potential well](). The escape rate is related to the ion temperature,
and both radial and axial escaping are considered.

### Charge exchange (CX)
The Charge exchange refers that the ion decreases the charge state because if the interaction with the neutral gases, could be described by:

$$B+A^{q+} = B^{n+}+A^{q-n+},$$

In the project, only the single and double charge exchange are considered, 
and the [Mueller-Salzborn formula](https://doi.org/10.1016/0375-9601(77)90672-7 "The charge exhcange cross section calculation") are used to estimate the cross section

## The ion temperature
In the following chapter, I will introduce the physcial progress that determines the ion temperature.

### The electron heating
The electron heating (Spitzer heating) comes from the collisional energy exchange between the electrons and ions:

$$\frac{dE_i}{dt}=\frac{4}{3}\nu_{ei}N_i\frac{m_e}{M_i}E_e,$$

The $\nu_{ei}=4\pi\frac{N_e}{v^3_e}(\frac{q_ie^2}{4\pi\epsilon_0m})ln\Lambda_i$ is the collision rate, where $v_e$ is the electron velocity and $ln\Lambda_i$ is the electron-ion Coulomb logarithm for the $i$ ion. The [electron-ion Coulomb logarithm](https://www.nrl.navy.mil/Portals/38/PDF%20Files/NRL_Formulary_2019.pdf?ver=p9F4Uq9wAtB0MPBwKYL9lw%3d%3d "Page34") is calculated based on $\frac{T_im_e}{m_i}$ criteria.

### The ion-ion energy exchange
The ion-ion energy exchanging progress is the key that determines the ion temperature. The energy exchange rate could be expressed as:

$$\sum_j\left(\frac{dE_i}{dt}\right)_j=\sum_j2\nu_{ij}N_i\frac{M_i}{M_j}\frac{kT_j-kT_i}{(1+\frac{M_ikT_j}{M_jkT_i})^{3/2}},$$

In the trap, the relativistic effects for the ion are not significant，so the collisional rate of the ion-ion are taken 
the [mixed ion-ion collisions situation](https://www.nrl.navy.mil/Portals/38/PDF%20Files/NRL_Formulary_2019.pdf?ver=p9F4Uq9wAtB0MPBwKYL9lw%3d%3d "Page34"):
$$\nu_{ii^\prime}=23-ln[\frac{ZZ^\prime(\mu+\mu^\prime)}{\mu T_{i^\prime}+\mu^{\prime}T_i}(\frac{n_iZ^2}{T_i}+\frac{n_{i^\prime}Z^{\prime 2}}{T_{i^\prime}})^{1/2}]$$

### The ion escape 
Both the contribution to ion temprature and the ions lossing rate are discussed in this chapter

The escaping ions also contribute the cooling of the ions. Combined the axial and radial escaping rate, the energy loss rate (solution of the Fokker-Plank equation) could be expressed as： 

$$\frac{dE_i^{ESC}}{dt}=-(\frac{2}{3}N_i\nu_i e^{-\omega_i}-R_i^{ESC})kT_i,$$

The $\nu_i$ is just the summation of the $\nu_{ij}$ collisional rate matrix between the ion and ion.

The escaped rate is taken:

$$\frac{dN_i}{dt}=\frac{-3}{\sqrt 2}N_i\nu_i\frac{e^{-\omega_i}}{\omega_i},$$

where $\omega_i=\frac{qV}{kT_i}$. For the radial situation, the magnetic field effect should be taken into consideration: $\omega_i=\frac{qV_r+2qBr_{DT}v_{\theta}}{kT_i}$, $v_{\theta}$ is the Maxwell averaged velocity of the ion.

## Code structure
In this chapter, I am going to introduce what you need to input to simulate the charge state distribution and the function of the methods. All the input parameters is defined at the properties part.

### ElectronBeamClass
This class is going to calculate the parameters that relates to the electron beam in EBIT.

#### What you need to Input
***Energy:*** The energy of the electron beam in electron volt (eV) unit;

***B_field:*** The axial magnetic field at the drift tubes part in tesla (T) unit;

***I_current:*** The emission current of the cathode in ampere (A) unit;

***Tc:*** The temperature of the cathode in keivn (K) unit;

***rc:*** The radius of the cathode emission part (crystal part) in metre (m) unit;

***Bc:*** The magentic field at the cathode in tesla (T) unit.

The above parameters are predefined, and you need to input it in the class before calling any method.
#### What is calculated in the properties
***rH:*** the calculated electron beam radius using herrman radius formula: 

$$r_H = r_B\sqrt{0.5+\sqrt{0.25+\frac{8m_ekt_cr_c^2}{eBr_B^2}+\frac{B_c^2r_c^4}{B^2r_B^4}}},$$

where the $r_B$ is the Brillouin Radius calculated by:

$$r_B=\sqrt{\frac{2m_eI_{emission}}{\pi\epsilon_0ev_eB^2}}.$$

***Velocity:*** the velocity of the electron in metre per second unit

***je:*** the current density at the trap center in ampere per square of centimeter ($A/cm^2$) unit

***jc:*** the current density at the cathode head in ampere per square of centimeter ($A/cm^2$) unit

***ne:*** the electron density at the trap center in per cube of centimetre ($cm^{-3}$)

#### Methods
***InitialParam(obj):*** The main function to calculate also the properties of the ElectronBeamClass. You need to call this function after creating the ElectronBeamClass and setting the trap parameters to the desired one. The following functions are both the sub-function of this one.

***rel_v():*** calculating the velocity of the electron beam considering the relativistic effect;

return type: ***scalar***

***BrillouinRadius():*** calculating the electron beam radius at the trap center using Brillouin formula;

return type: ***scalar***

***HerrmanRadius():*** calculating the electron beam radius at the trap ceneter using Herrman formula；

return type: ***scalar***

***ElectronDensity():*** calcuating the area density and the density of the electron beam at the trap center；

return type: ***scalar***

***ElectronDensityCathode():*** calculating the area density of the electron beam at the cathode head.

return type: ***scalar***

### ElementClass
This class is going to calculate all kinds of cross sections(RR, EII, CX, DCX) for a specific element.

#### What you need to input
When simulating the charge state distribution and the ion temperature in EBIT, you typically need to create two elements class: one is the target element, another is the light cooling element(evaporative cooling: baground gas also takes part in). In the following context, these two will be called HCI and LCI.

***Zl:*** the proton number of the LCI in atomic unit (amu);

***Ml:*** the mass number of the LCI in atomic unit (amu);

***Path_l:*** the ionization potential file for the LCI (in keV unit in the file);

***Path:*** the ionization potential file for the HCI (in keV unit in the file);

***Zh:*** the proton number of the HCI in atomic unit (amu);

***Mh:*** the mass number of the HCI in atomic unit (amu);

***Path:*** the ionization potential file for the HCI (in keV unit in the file);

#### What is calculated in the properties
***EII:*** the electron impact ionization cross section for the specific element in $cm^{2}$ unit;

***RR:*** the radiative recombination cross section for the specific element in $cm^{2}$ unit;

***CX:*** the charge exchange cross section for the specific element in $cm^{2}$ unit;

***DCX:*** the double electron charge exchange for the specific element in $cm^{2}$ unit;

***ion_pot:*** the ionization potential array for each HCI ionic species in eV uint;

***l_ion_pot:*** the ionization potential array for each LCI ionic species in eV uint;

***eBeam:*** the reference ElectronBeamClass 

#### Methods
***InitialPara(obj,Z,mh):*** the initialization of the ElementClass, here the Z and mh are the atomic number and mass number for HCI, respectively. You also need to call this function in the main script after the paramter setting.

***ion_potload(Path):*** the loading and unit changing of ionization potential for each ionic species.

return type ***(Z+1)x2 vector***

***EII_cross():*** the electron ionization cross section calculation in $cm^2$ unit.

return type ***(Z+1)x1 vector***

***RR_cross():*** the radiative recombination cross section calculation in $cm^2$ unit.

return type ***(Z+1)x1 vector***

***CX_cross, DCX_cross():*** the single and double electron charge exchange cross section calculaiton in $cm^2$ unit

return type ***(Z+1)x1 vector***

### PlasmaClass
This class will construct all kinds of plasma parameters and two-version equation rate equation (with and without ion temperature calculation)

#### What you need to input
***r_DT:*** The center drift tube radius in meter (m) unit;

***Vw:*** The potential applied on the surrending drift tubes in Volt (V) unit (The approximation is based on the square potential well);

***Ti_min:*** The minimum ion temperature in electron volt (eV). (The temperature can not be negative, so a minimum ion temperature is set)

***ni_min:*** The minimum ion density that taking part in the plasma evolution, once the density of the specific ion species is lower than ni_min, the changing rate the ion will be set to 0 ($cm^{-3}$);

***n0_HCI***: The absolute initial density of the neutral HCI atoms density ($cm^{-3}$);

***n0_LCI***: The absolute initial density of the neutral LCI atoms density ($cm^{-3}$);

***S_term***: Source term of the continous injected element. Please use it in the no-temperature situation. If you need to use it in ion evolution with temperature, please set the n0_HCI = 0, and do not set the S_term value too large (1e+4 is a suitable value).

#### What is calculated in this class
***f:*** The overlapping factor beween the electron and the ion cloud. In the current situation, *f* is setting to 1 for each ionic species;

***Phi_r:*** The radial trap potential from the electron beam. Now, the trap potential calculaiton is based on the cylindrical approximation of the electron beam and the Non-relativistic calculation of the line charge density (V);

***HCI:*** The reference class of the heavy target HCI ElementClass;

***LCI:*** The reference class of the light cooling LCI ElementClass;

***Ti:*** The ion temperature vector for all HCI and LCI (eV);

***dTidt:*** The first time derivative of the ion temperature (eV/s);

***Te:*** The electron temperature (eV). The electron temperature calculation is approximated as $0.1\frac{j_e}{j_c}$;

***ni:*** The ion density ($cm^{-3}$);

***dnidt:*** The first time derivative of the ion density ($cm^{-3}/s$);

***eBeam:*** The reference class of the ElectronBeamClass;

***ei_e_rate:*** The electron ion energy exchange rate ($eV.cm^3/s$), which is the $\nu_{ei}$;

***clog_ij:*** The ion-ion coulomb logarithum martix;

***nu_ij:*** The ion-ion collision rate between the i and j ionic species;

***ii_ex_rate:**** The ion-ion energy exchange rate between the i ionic specie and all other ions;

***omega_i:*** The characteristic frequency of the ion escaping from the axial direction only taking axial potential into consideration;

***omega_r_i:*** The characteristic frequency of the ion escaping from the radial direction combined the Beam potential with the magentic field;

***esc_e_rate:*** The 
