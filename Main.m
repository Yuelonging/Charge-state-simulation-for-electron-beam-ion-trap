%%
clc
clear
%%
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

%%
Y0 = transpose([EBIT.ni,EBIT.Ti]);
t_span = [0,100];
opts = odeset('NonNegative', [1:2*(Zh+Zl+2)], 'MaxStep',2);
ode_fun = @(t, Y) ode_wrapper(t, Y, EBIT);
[t, Y] = ode15s(ode_fun,t_span, Y0,opts);
%%
figure(1)
semilogx(t,Y(:,Zh-5:Zh+1)./sum(Y(:,1:Zh+1),2),'LineWidth',1.5);
ylim([0,0.6]);
xlim([1e-1,1e+2]);
grid on
legend('C','B','Be','Li','He','H','Bare');
title("^{120}_{50}Sn Charge state E_e = 50000 eV ")
xlabel("Time (s)")
%%
figure(2)
id = Zh+Zl+2;
semilogx(t,Y(:,(1:Zh+1)+id),'.');
% ylim([1e+2,1e+10]);
xlim([1e-3,1e+3]);

%%
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