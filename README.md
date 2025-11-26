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

***ii_ex_rate:*** The ion-ion energy exchange rate between the i ionic specie and all other ions;

***omega_i:*** The characteristic frequency of the ion escaping from the axial direction only taking axial potential into consideration;

***omega_r_i:*** The characteristic frequency of the ion escaping from the radial direction combined the Beam potential with the magentic field;

***esc_e_rate:*** The energy loss rate in the axial direction ($eV.cm^3/s$);

***esc_e_r_rate:*** The energy loss rate in the radial direction ($eV.cm^3/s$);

***RR_rate:*** The RR rate, which is $\frac{j_e}{e}n_i\sigma_{RR}f$;

***CX_rate:*** The single electron charge exchange rate, which is $n_0n_i\sigma_{CX}\overline{v_i}$. $\overline{v_i}=\sqrt{\frac{8T_ie}{m_i\pi}}$;

***DCX_rate:*** The double electron charge exchange rate. Formula same as the CX_rate;

***EII_rate:*** The EII rate, which is $\frac{j_e}{e}n_i\sigma_{EII}f$;

***ESC_a_rate:*** The escape rate of the ion from the trap axially;

***ESC_r_rate:*** The escape rate of the ion from the trap radially.

#### Methods
***InitialParm(obj):*** Initilization of the Te, Ti, dnidt, dTidt, f, Phi_r.

&emsp;***Te_calculation():*** Electron temperature estimation.

&emsp;return type ***scalar***

&emsp;***Ti_initialization():*** Initialization of the ion temperature and temperature changing rate. The Ti is set to 1.1*Ti_min, and the dTidt is set to 0.

&emsp;return type Ti:***1x(Zh+Zl+2) Vector*** dTidt:***1x(Zh+Zl+2) Vector***

&emsp;***ni_initialization():*** Initialization of the ion ambulance, the HCI I and LCI I is set to a given value.

&emsp;return type ni:***1x(Zh+Zl+2) Vector*** dnidt:***1x(Zh+Zl+2) Vector***

&emsp;***f_Initialization()*** Initialization of the overlapping factor for the electron beam and ion beam.

&emsp;return type: ***1x(Zh+Zl+2)***

&emsp;***radial_trap_strength()*** Directly calculate the potential difference between the center beam and trap edge.

&emsp;return type: ***scalar***

***c_ei():*** Calculation for the coulomb logarithum for electron and ion, and the energy exchange rate between ion and electron.

return type: ei_e_rate:***1x(Zh+Zl+2) vector***

***ii_ex_r():*** Calculation for the ion-ion energy exchange rate.

return type: ii_ex_rate:***1x(Zh+Zl+2) vector***

&emsp;***clog_ii(Index)*** The index here is the a traversal for the [HCI, LCI] array (corresponding the charge state [I,II,III...N, I,II,III...M]).

&emsp;return type: ***1x(Zh+Zl+2) vector***, corresponding i ion specie between all the other ion.

&emsp;***nu_ii(clog_ij_r):*** clog_ij_r: the coulomb logarithum matrix after the index loop, so it should be a ***(Zh+Zl+2)x(Zh+Zl+2) matrix***. The function

&emsp;will calculate the ion-ion collisional rate.

&emsp;return type:***(Zh+Zl+2)x(Zh+Zl+2) matrix***

&emsp;***ii_ex():*** transpose and sum the nu_ij matrix, and calculate the energy transfer rate.

&emsp;return type:***1x(Zh+Zl+2) vector***

***esc_e_r(obj), esc_e_r_rate(obj):*** The energy loss rate from the trap radially and axially. The radial lossing rate including the magentic field contribution.

return type: ***1x(Zh+Zl+2) vector***

***RR_r(obj),EII_r(obj),CX_r(obj),DCX_r(obj),ESC_a_r(obj),ESC_r_r(obj):*** The calculation of the ion RR EII CX DCX ESC rate ($/s/cm^3$)

return type: ***1x(Zh+Zl+2) vector***

***Rate_eq(obj):*** Construct the charge state and temperature evolution equations

$$
\begin{align}
\frac{dn_0}{dt}=S_{term}+\frac{j_e}{e}f(-n_0\sigma^{EII}_{0&rarr;1}+n_1\sigma^{RR}_{1&rarr;0})+N_0n_1\overline{v_1}\sigma^{CX}_{1&rarr;0},\\
...\\
\frac{dn_i}{dt}=\frac{j_e}{e}f(-n_i\sigma^{EII}_{i&rarr;i+1}+n_{i-1}\sigma^{EII}_{i-1&rarr;i}+n_{i+1}\sigma^{RR}_{i+1&rarr;i}-n_i\sigma^{RR}_{i&rarr;i-1})
+N_0(n_{i+1}\overline{v}_{i+1}\sigma^{CX}_{i+1&rarr;i}-n_{i}\overline{v}_{i}\sigma^{CX}_{i&rarr;i-1})
-\frac{3}{\sqrt{2}}n_i\nu_i\frac{e^{-\omega_i}}{\omega_i},\\
...\\
\frac{dn_m}{dt}=\frac{j_e}{e}f(n_{m-1}\sigma^{EII}_{m-1&rarr;m}-n_m\sigma^{RR}_{m&rarr;m-1})
-N_0(n_{m}\overline{v}_{m}\sigma^{CX}_{m&rarr;m-1})
-\frac{3}{\sqrt{2}}n_m\nu_m\frac{e^{-\omega_m}}{\omega_m},\\
\frac{dT_i}{dt}=\frac{2}{3}n_e\frac{m_e^2}{m_iv_e}4\pi(\frac{q_ie^2}{(4\pi\epsilon_0m_e})^2ln\Lambda_i
+2\nu_{ij}\frac{M_i}{M_j}\frac{T_j-T_i}{(1+\frac{M_iT_j}{M_jT_i})^{3/2}}
-\frac{2}{3}\nu_ie^{-\omega_i}T_i
\end{align}
$$

***Rate_eq_without_temperature(obj):*** To construct the rate equation without the ion temperature parameter taking part in the charge state ecolution

$$
\begin{align}
\frac{dn_0}{dt}=S_{term}+\frac{j_e}{e}f(-n_0\sigma^{EII}_{0&rarr;1}+n_1\sigma^{RR}_{1&rarr;0})+N_0n_1\overline{v_1}\sigma^{CX}_{1&rarr;0},\\
...\\
\frac{dn_i}{dt}=\frac{j_e}{e}f(-n_i\sigma^{EII}_{i&rarr;i+1}+n_{i-1}\sigma^{EII}_{i-1&rarr;i}+n_{i+1}\sigma^{RR}_{i+1&rarr;i}-n_i\sigma^{RR}_{i&rarr;i-1})
+N_0(n_{i+1}\overline{v}_{i+1}\sigma^{CX}_{i+1&rarr;i}-n_{i}\overline{v}_{i}\sigma^{CX}_{i&rarr;i-1})\\
...\\
\frac{dn_m}{dt}=\frac{j_e}{e}f(n_{m-1}\sigma^{EII}_{m-1&rarr;m}-n_m\sigma^{RR}_{m&rarr;m-1})
-N_0(n_{m}\overline{v}_{m}\sigma^{CX}_{m&rarr;m-1})
\end{align}
$$

### Main
In this chapter, I will introduce how to construct a script. The progress is based on that you have set all the required parameters in the class (There is an instruction in each chapter ***What you need to input*** section).
```
Ebeam = ElectronBeamClass();
Ebeam = Ebeam.InitialParam();
Pb = ElementClass();
Pb.eBeam = Ebeam;
Pb.Path = "C:\Users\jialin\Desktop\phdWorks\EBIT\Script\ChargeSimulation\U.ion";
Pb = Pb.InitialParam(92,238);
Ne = ElementClass();
Ne.eBeam = Ebeam;
Ne.Path = "C:\Users\jialin\Desktop\phdWorks\EBIT\Script\ChargeSimulation\Ne.ion";
Ne = Ne.InitialParam(10,20);
```
First, you need to the ElectronBeamClass, ElementClass, loading the Ionization potential for both HCI and LCI, and call the ***InitialParam()*** to calculate the needed property in the class.
```
EBIT=PlasmaClass();
EBIT.HCI=Pb;EBIT.LCI=Ne;EBIT.eBeam=Ebeam;
EBIT = EBIT.InitialParm();
EBIT.ei_e_rate = EBIT.c_ei();
EBIT.ii_ex_rate=EBIT.ii_ex_r();
EBIT.esc_e_rate=EBIT.esc_e_r();
EBIT.esc_e_r_rate=EBIT.esc_e_r_r();
EBIT.RR_rate=EBIT.RR_r();
EBIT.EII_rate=EBIT.EII_r();
EBIT.CX_rate=EBIT.CX_r();
EBIT.DCX_rate=EBIT.DCX_r();
EBIT.ESC_a_rate=EBIT.ESC_a_r();
EBIT.ESC_r_rate=EBIT.ESC_r_r();
[EBIT.dnidt,EBIT.dTidt] = EBIT.Rate_eq();
Zh = EBIT.HCI.Zh;
Zl = EBIT.HCI.Zl;
```
Calling the PlasmaClass and intializae all the rate coeficient.
```
Y0 = transpose([EBIT.ni,EBIT.Ti]);
t_span = [0,100];
opts = odeset('NonNegative', [1:2*(Zh+Zl+2)], 'MaxStep',2);
ode_fun = @(t, Y) ode_wrapper(t, Y, EBIT);
[t, Y] = ode23t(ode_fun,t_span, Y0,opts);
```
Set the odeset, and solve the equation
```
function dYdt = ode_wrapper(t, Y, obj)
    Zh = obj.HCI.Zh;
    Zl = obj.HCI.Zl;
    obj.ni = transpose(Y(1:Zh+Zl+2));
    obj.Ti = transpose(Y(Zh+Zl+3:end));
    eq_updating(obj);
    dYdt = transpose([obj.dnidt,obj.dTidt]);
    function eq_updating(obj)
        obj.ni(find(obj.ni<obj.ni_min))=obj.ni_min;
        obj.Ti(find(obj.Ti<obj.Ti_min))=obj.Ti_min;
        obj.ei_e_rate = obj.c_ei();
        obj.ii_ex_rate=obj.ii_ex_r();
        obj.esc_e_rate=obj.esc_e_r();
        obj.esc_e_r_rate=obj.esc_e_r_r();
        obj.RR_rate=obj.RR_r();
        obj.EII_rate=obj.EII_r();
        obj.CX_rate=obj.CX_r();
        obj.DCX_rate=obj.DCX_r();
        obj.ESC_a_rate=obj.ESC_a_r();
        obj.ESC_r_rate=obj.ESC_r_r();
        [obj.dnidt,obj.dTidt] = obj.Rate_eq();
    end
end
```
The ode wrapper, just updating all the rate coefficient after a time step.

### Main_no_temperature
Here is an version without temperature calculation.
```
Ebeam = ElectronBeamClass();
Ebeam = Ebeam.InitialParam();
Pb = ElementClass();
Pb.eBeam = Ebeam;
Pb.Path = "C:\Users\jialin\Desktop\phdWorks\EBIT\Script\ChargeSimulation\U.ion";
Pb = Pb.InitialParam(92,238);
Ne = ElementClass();
Ne.eBeam = Ebeam;
Ne.Path = "C:\Users\jialin\Desktop\phdWorks\EBIT\Script\ChargeSimulation\Ne.ion";
Ne = Ne.InitialParam(10,20);
%% Initialization of the rate equation
EBIT=PlasmaClass();
EBIT.HCI=Pb;EBIT.LCI=Ne;EBIT.eBeam=Ebeam;

%%
EBIT = EBIT.InitialParm();
EBIT.RR_rate=EBIT.RR_r();
EBIT.EII_rate=EBIT.EII_r();
EBIT.CX_rate=EBIT.CX_r();

[EBIT.dnidt,EBIT.dTidt] = EBIT.Rate_eq_without_temperature();
Zh = EBIT.HCI.Zh;
Zl = EBIT.HCI.Zl;

%%
Y0 = transpose([EBIT.ni,EBIT.Ti]);
t_span = [0,1000];
opts = odeset('NonNegative', [1:2*(Zh+Zl+2)], 'MaxStep',2);
ode_fun = @(t, Y) ode_wrapper(t, Y, EBIT);
[t, Y] = ode23t(ode_fun,t_span, Y0,opts);

function dYdt = ode_wrapper(t, Y, obj)
    Zh = obj.HCI.Zh;
    Zl = obj.HCI.Zl;
    obj.ni = transpose(Y(1:Zh+Zl+2));
    obj.Ti = transpose(Y(Zh+Zl+3:end));
    eq_updating(obj);
    dYdt = transpose([obj.dnidt,obj.dTidt]);
    function eq_updating(obj)
        obj.ni(find(obj.ni<obj.ni_min))=obj.ni_min;
        obj.Ti(find(obj.Ti<obj.Ti_min))=obj.Ti_min;
        obj.RR_rate=obj.RR_r();
        obj.EII_rate=obj.EII_r();
        obj.CX_rate=obj.CX_r();
        [obj.dnidt,obj.dTidt] = obj.Rate_eq_without_temperature();
    end
end
```
The Initialization progress is the same. Just call the rate_eq_without_temperature function instead of rate_eq. And there is a modification for the wrapper function.
