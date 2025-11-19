classdef PlasmaClass < handle
    properties
        %Parameters:
        %f: The overlap factor between the ion and trap (maxium 1)
        %Phi_r (V) The radial trap potential from the electron beam
        %r_DT (m) The center drift tube radius
        %Vw (V) the potential well applied on drift tubes of to trap the ion
        %HCI the highly charged ions EI,RR,CX
        %LCI the ions with a relative low chagred ions, participate in the
        %cooling and CX(Ne I, for the Mueller Salzborn formula)
        %Ti (eV) Ion temperaure
        %dTidt (eV/s) The first time derivative of the ion temperature
        %Te (eV) electron temeperatue
        %Ti_min (eV) The allowed minimun ion temperature
        %ni (/cm^3) The ion density for HCI and LCI
        %dnidt (/cm^3/s) The first time derivative of the ion density 
        %ni_min (/cm^3) The allowed minimun Ion density
        %n0_HCI (/cm^3) The initial neutral density of the HCI atoms
        %n0_LCI (/cm^3) The initial neutral density of the LCI atoms
        %eBeam: the electron beam struct
        %ei_e_rate (eV.cm^3/s) The electron ion energy exchange rate
        %clog_ei: The electron-ion coulomb logarithum
        %clog_ij: The ion-ion coulomb logarthum
        %nu_ij: The ion-ion collision rate between the i and j ionic species
        %ii_ex_rate (eV.cm^3/s) The ion ion exchange rate between the i and
        %all other ions
        %omega_i: The potential over Ti characteristic frequency
        %omega_r_i: Same as omega_i but for radial
        %esc_e_rate (eV.cm^3/s) The energy loss rate from the trap
        f;
        Phi_r;
        r_DT = 2e-3;
        Vw = 100;
        HCI;
        LCI;
        Ti;
        dTidt;
        Te;
        Ti_min = 0.026;
        ni;
        dnidt
        ni_min = 1e-6;
        n0_HCI = 1e+8;
        n0_LCI = 1e+6;
        eBeam;
        ei_e_rate;
        clog_ij;
        nu_ij;
        ii_ex_rate;
        omega_i;
        omega_r_i;
        esc_e_rate;
        esc_e_r_rate;
        RR_rate;
        CX_rate;
        DCX_rate;
        EII_rate;
        ESC_a_rate;
        ESC_r_rate;
    end
    methods
        function obj=InitialParm(obj)
            obj.Te = Te_calculation();
            [obj.Ti,obj.dTidt] = Ti_Initialization();
            [obj.ni,obj.dnidt] = ni_Initialization();
            obj.f = f_Initalizaion();
            obj.Phi_r = radial_trap_strength();

            function Te = Te_calculation()
                %The Calculation of the electron temperatue
                %ref: https://indico.cern.ch/event/199158/contributions/
                % 379787/attachments/295625/413143/T_II_2_Marrs.pdf
                Te = 0.1*obj.eBeam.je/obj.eBeam.jc;
            end
    
            function [Ti,dTidt] = Ti_Initialization()
                %All the ion temprature is set to 300 K
                Ti = ones(1,obj.HCI.Zh+obj.LCI.Zh+2)*obj.Ti_min*1.1;
                dTidt=zeros(1,obj.HCI.Zh+obj.LCI.Zh+2)*obj.Ti_min*1.1;
            end

            function [ni,dnidt] = ni_Initialization()
                %All the ion are set as The 
                ni = ones(1,obj.HCI.Zh+obj.LCI.Zh+2)*obj.ni_min*1.1;
                ni(1) = obj.n0_HCI;
                ni(obj.HCI.Zh+2) = obj.n0_LCI;
                dnidt = zeros(1,obj.HCI.Zh+obj.LCI.Zh+2)*obj.ni_min;
            end

            function f = f_Initalizaion()
                f= ones(1,obj.HCI.Zh+obj.LCI.Zh+2);
            end

            function Phi_r = radial_trap_strength()
                %The calculation of the potential difference between the
                %center beam and the trap edge
                rtrap = obj.r_DT;
                rbeam = obj.eBeam.rH;
                ve = obj.eBeam.Velocity;
                eps0 = 8.854e-12;
                Ie = obj.eBeam.I_current;
                Phi_r = Ie/(4*pi*eps0*ve)*(2*log(rtrap/rbeam)+1);
            end

        end

        function ei_e_rate=c_ei(obj)
            %Here is a calculation for the coulomb logarthium for electron
            %ion and the energy exchange rate between ion and electron
            %Page 34
            %ref: https://www.nrl.navy.mil/Portals/38/PDF%20Files/
            % NRL_Formulary_2019.pdf?ver=p9F4Uq9wAtB0MPBwKYL9lw%3d%3d
            %equation (7)
            %ref: Phys. Rev. A., 43. 4861
            coulog_ei=[];
            ne = obj.eBeam.ne;
            me = obj.eBeam.me_kg;
            mi_H = obj.HCI.Mh*1.6605e-27;
            mi_L = obj.LCI.Mh*1.6605e-27;
            Zh = obj.HCI.Zh; Zl =obj.HCI.Zl;

            %The calculation of the HCI part
            for  i = 1:Zh+1
                q=i-1;%The neutral atom with charge state 0
                if obj.Ti(i)*me/mi_H<obj.Te&&obj.Te<10*q^2
                    coulog_ei(i) = 23-log(ne^0.5*q*obj.Te^(-1.5));
                elseif obj.Ti(i)*me/mi_H<10*q^2&&10*q^2<obj.Te
                    coulog_ei(i) = 24-log(ne^0.5/obj.Te);
                elseif obj.Te<obj.Ti(i)*me/mi_H
                    coulog_ei(i) = 16-log(obj.ni(i)^0.5*obj.Ti(i)^(-1.5)*q^2*obj.HCI.Mh);
                end
            end
            %The calculation of the LCI part
            for i = Zh+2:Zh+Zl+2
                q=i-Zh-2;
                 if obj.Ti(i)*me/mi_L<obj.Te&&obj.Te<10*q^2
                    coulog_ei(i) = 23-log(ne^0.5*q*obj.Te^(-1.5));
                elseif obj.Ti(i)*me/mi_L<10*q^2&&10*q^2<obj.Te
                    coulog_ei(i) = 24-log(ne^0.5/obj.Te);
                elseif obj.Te<obj.Ti(i)*me/mi_L
                    coulog_ei(i) = 16-log(obj.ni(i)^0.5*obj.Ti(i)^(-1.5)*q^2*obj.LCI.Mh);
                end
            end
            
            %Calculate the coulomb logauthrim cross section for each charge
            %sate
            q = [0:Zh, 0:Zl];
            e = 1.602e-19;
            eps0 = 8.854e-12;
            % Y = 4*pi*(q*e^2/(4*pi*eps0*me)).^2.*coulog_ei;
            Y = 4*pi*(q*e^2/(4*pi*eps0*me)).^2.*coulog_ei/obj.eBeam.Velocity^4;
            Mi = [ones(1,Zh+1)*mi_H, ones(1,Zl+1)*mi_L];
            % ei_e_rate = 2/3*Y*ne*1e+6*me^2./Mi/(obj.eBeam.Velocity);
            % ei_e_rate = ei_e_rate/e;%Because of the kbTi is in (eV)
            ei_e_rate = 2/3*Y*ne*1e+6*obj.eBeam.Velocity*2*me./Mi*obj.eBeam.Energy;
            C1 = (obj.ni>=obj.ni_min);
            ei_e_rate = ei_e_rate.*C1;
        end

        function ii_ex_rate=ii_ex_r(obj)
            %The calculation of the ion-ion energy exchange rate
            clog_ij_r=[];
            Index = 1:obj.HCI.Zh+obj.LCI.Zh+2;
            for i = 1:length(Index)
                clog_ij_r = [clog_ij_r;clog_ii(Index(i))];
            end
            clog_ij_r(find(clog_ij_r==Inf))=0;
            obj.clog_ij = clog_ij_r;
            obj.nu_ij=nu_ii(clog_ij_r);
            ii_ex_rate = ii_ex();


            function clog_ij=clog_ii(Index)
                %The calculation for the coulomb logarithm between ion and
                %ion
                %Index; the Index to find the corrsponding ion density
                %temperature
                %Page 34
                %ref: https://www.nrl.navy.mil/Portals/38/PDF%20Files/
                % NRL_Formulary_2019.pdf?ver=p9F4Uq9wAtB0MPBwKYL9lw%3d%3d
                if Index<=obj.HCI.Zh+1
                    %The charge Index need to plus 1 to match the Index
                    %1->Ne I; 2->Ne II
                    q=Index-1;
                    q_prime = [0:obj.HCI.Zh,0:obj.LCI.Zh];
                    M = obj.HCI.Mh; M_prime = obj.HCI.Mh;
                    %The coulomb collision with the HCI
                    clog_ij_H=23-log(q*q_prime*(M+M_prime)./(M*obj.Ti+M_prime*obj.Ti(q+1))...
                        .*sqrt(obj.ni(Index)*q.^2/obj.Ti(Index)+obj.ni.*(q_prime).^2./obj.Ti));
                    %The coulomb collision with the LCI
                    M = obj.HCI.Mh; M_prime = obj.HCI.Ml;
                    clog_ij_L=23-log(q*q_prime*(M+M_prime)./(M*obj.Ti+M_prime*obj.Ti(q+1))...
                        .*sqrt(obj.ni(Index)*q.^2/obj.Ti(Index)+obj.ni.*(q_prime).^2/obj.Ti));
                    clog_ij = [clog_ij_H(1:obj.HCI.Zh+1),clog_ij_L(end-obj.LCI.Zh:end)];
                elseif Index<=obj.HCI.Zh+obj.LCI.Zl+2&&Index>obj.HCI.Zh+1
                    q=Index-obj.HCI.Zh-2;
                    q_prime = [0:obj.HCI.Zh,0:obj.LCI.Zh];
                    M = obj.HCI.Mh; M_prime = obj.HCI.Mh;
                    %The coulomb collision with the HCI
                    clog_ij_H=23-log(q*q_prime*(M+M_prime)./(M*obj.Ti+M_prime*obj.Ti(q+1))...
                        .*sqrt(obj.ni(Index)*q.^2/obj.Ti(Index)+obj.ni.*(q_prime).^2/obj.Ti));
                    %The coulomb collision with the LCI
                    M = obj.HCI.Mh; M_prime = obj.HCI.Ml;
                    clog_ij_L=23-log(q*q_prime*(M+M_prime)./(M*obj.Ti+M_prime*obj.Ti(q+1))...
                        .*sqrt(obj.ni(Index)*q.^2/obj.Ti(Index)+obj.ni.*(q_prime).^2/obj.Ti));
                    clog_ij = [clog_ij_H(1:obj.HCI.Zh+1),clog_ij_L(end-obj.LCI.Zh:end)];
                end
            end

            function nu_ij=nu_ii(clog_ij_r)
                %The calculation for the ion ion collisional rate
                %clog_ij the ion-ion coulomb logarithm matrix
                %euqation (9)
                %ref: Phys. Rev. A., 43. 4861
                Nj = obj.ni;
                Zi = transpose([0:obj.HCI.Zh, 0:obj.LCI.Zh]);
                Zj = [0:obj.HCI.Zh, 0:obj.LCI.Zh];
                e = 1.6022e-19;
                eps0 = 8.854e-12;
                Mi = transpose([ones(1,obj.HCI.Zh+1)*obj.HCI.Mh*1.6605e-27,...
                    ones(1,obj.LCI.Zh+1)*obj.LCI.Mh*1.6605e-27]);%kg
                %Here the Nj times 1e+6 to ensure the SI unit
                nu_ij = 4/3*sqrt(2*pi)/(4*pi*eps0)^2.*(Nj*1e+6).*(Zi.*Zj*e^2./Mi).^2.*(Mi./(transpose(obj.Ti)*e)).^1.5.*clog_ij_r;
                %Ensure that the too low and too samll energy ion do not join the exchange 
                C1 = transpose(obj.ni>=obj.ni_min)*(Nj>=obj.ni_min);
                C2 = transpose(obj.Ti>=obj.Ti_min)*(obj.Ti>=obj.ni_min);
                nu_ij = nu_ij.*C1.*C2;                
            end

            function ii_ex_rate=ii_ex()
                %The calculation for the ion ion energy exchange equation
                Mi = transpose([ones(1,obj.HCI.Zh+1)*obj.HCI.Mh,...
                    ones(1,obj.LCI.Zh+1)*obj.LCI.Mh]);
                Mj = [ones(1,obj.HCI.Zh+1)*obj.HCI.Mh,...
                    ones(1,obj.LCI.Zh+1)*obj.LCI.Mh];
                Tj = obj.Ti;
                ii_ex_rate = 2*obj.nu_ij.*Mi./Mj.*(Tj-transpose(obj.Ti))./(1+Mi.*Tj./(Mj.*transpose(obj.Ti))).^1.5;
                %for each j, we sum the coloum
                ii_ex_rate = transpose(sum(ii_ex_rate,2));
            end
        end

        function esc_e_rate = esc_e_r(obj)
            %The calculation of the energy loss rate from the trap
            Zi = [0:obj.HCI.Zh, 0:obj.LCI.Zh];
            obj.omega_i = Zi*obj.Vw./obj.Ti;%Here Ti is in eV, so there is no e on the top;
            nu_i = transpose(sum(obj.nu_ij, 2));
            esc_e_rate = -(2/3).*nu_i.*exp(-obj.omega_i).*obj.Ti;
            C1 = obj.ni>obj.ni_min;
            esc_e_rate = esc_e_rate.*C1;
        end

        function esc_e_r_rate = esc_e_r_r(obj)
            %The calculation of the energy loss rate from the trap
            au2kg = 1.6605e-27;
            e = 1.602e-19;
            Zi = [0:obj.HCI.Zh, 0:obj.LCI.Zh];
            Zh = obj.HCI.Zh; Zl = obj.LCI.Zh;
            Mh = obj.HCI.Mh; Ml = obj.LCI.Mh;
            B = obj.eBeam.B_field;
            Mi = [ones(1,Zh+1)*Mh*au2kg,ones(1,Zl+1)*Ml*au2kg];%kg
            v_theta = sqrt(2*obj.Ti*e./(3*Mi));
            obj.omega_r_i = (Zi*obj.Phi_r+2*Zi*B*obj.r_DT.*v_theta)./(obj.Ti);
            nu_i = transpose(sum(obj.nu_ij, 2));
            esc_e_r_rate = -(2/3).*nu_i.*exp(-obj.omega_r_i).*obj.Ti;
            C1 = obj.ni>obj.ni_min;
            esc_e_r_rate = esc_e_r_rate.*C1;
        end

        function RR_rate = RR_r(obj)
            sigma_RR = [obj.HCI.RR,obj.LCI.RR];
            e = 1.602e-19;
            je = obj.eBeam.je;
            RR_rate = je/e*obj.ni.*sigma_RR.*obj.f;
        end

        function EII_rate = EII_r(obj)
            sigma_EII = [obj.HCI.EII,obj.LCI.EII];
            e = 1.602e-19;
            je = obj.eBeam.je;
            EII_rate = je/e*obj.ni.*sigma_EII.*obj.f;
        end

        function CX_rate = CX_r(obj)
            sigma_CX = [obj.HCI.CX,obj.LCI.CX];
            e = 1.602e-19;
            au2kg = 1.6605e-27;
            Zh = obj.HCI.Zh;
            Zl = obj.HCI.Zl;
            Mh = obj.HCI.Mh;
            Ml = obj.LCI.Mh;
            vi_bar = vi_b();
            n0 = obj.ni(Zh+2);
            CX_rate = n0*obj.ni.*sigma_CX.*vi_bar;
            function vi_bar=vi_b()
                Mi = [ones(1,Zh+1)*Mh*au2kg,...
                    ones(1,Zl+1)*Ml*au2kg];%kg
                vi_bar = sqrt(8*obj.Ti*e./Mi/pi)*100;
                %change the m/s to cm/s to meet the unit
            end
        end

        function DCX_rate = DCX_r(obj)
            sigma_DCX = [obj.HCI.DCX,obj.LCI.DCX];
            e = 1.602e-19;
            au2kg = 1.6605e-27;
            Zh = obj.HCI.Zh;
            Zl = obj.HCI.Zl;
            Mh = obj.HCI.Mh;
            Ml = obj.LCI.Mh;
            vi_bar = vi_b();
            n0 = obj.ni(Zh+2);
            DCX_rate = n0*obj.ni.*sigma_DCX.*vi_bar*0;
            %Here we set DCX rate to 0 to close it
            function vi_bar=vi_b()
                Mi = [ones(1,Zh+1)*Mh*au2kg,...
                    ones(1,Zl+1)*Ml*au2kg];%kg
                vi_bar = sqrt(8*obj.Ti*e./Mi/pi)*100;
                %change the m/s to cm/s to meet the unit
            end
        end

        function ESC_a_rate = ESC_a_r(obj)
            nu_i = transpose(sum(obj.nu_ij,2));
            ESC_a_rate = -3/sqrt(2)*obj.ni.*nu_i.*exp(-obj.omega_i)./obj.omega_i;
            % C1 = obj.ni>obj.ni_min;
            % ESC_a_rate = C1.*ESC_a_rate;
            ESC_a_rate(1) = 0; ESC_a_rate(obj.HCI.Zh+2) = 0;
        end

        function ESC_r_rate = ESC_r_r(obj)
            nu_i = transpose(sum(obj.nu_ij,2));
            ESC_r_rate = -3/sqrt(2)*obj.ni.*nu_i.*exp(-obj.omega_r_i)./obj.omega_r_i;
            ESC_r_rate(1) = 0; ESC_r_rate(obj.HCI.Zh+2) = 0;
        end

        function [dnidt,dTidt]= Rate_eq(obj)
            dnidt = ones(1,obj.HCI.Zh+obj.HCI.Zl+2);
            dTidt = obj.ei_e_rate + obj.ii_ex_rate + obj.esc_e_rate+obj.esc_e_r_rate;
            Zh = obj.HCI.Zh;
            EII_im1_i_h = obj.EII_rate(1:Zh);
            EII_i_ip1_h = obj.EII_rate(2:(Zh+1));
            RR_ip1_i_h = [obj.RR_rate(3:(Zh+1)),0];
            RR_i_im1_h = obj.RR_rate(2:(Zh+1));
            CX_ip1_i_h = [obj.CX_rate(3:(Zh+1)),0];
            CX_i_im1_h = obj.CX_rate(2:(Zh+1));
            ESC_a_h = obj.ESC_a_rate(2:(Zh+1));
            ESC_r_h = obj.ESC_r_rate(2:(Zh+1));
            % dnidt(1) = -obj.EII_rate(1)+obj.RR_rate(2)+obj.CX_rate(2)+obj.DCX_rate(3); 
            dnidt(1) = 0;


            dnidt(2:(Zh+1)) = EII_im1_i_h-EII_i_ip1_h+RR_ip1_i_h-RR_i_im1_h+CX_ip1_i_h-CX_i_im1_h+ESC_a_h+ESC_r_h;
            DCX_ip2_i = [obj.DCX_rate(4:Zh+1),0,0];
            DCX_i_im2 = [0,obj.DCX_rate(3:Zh+1)];

            dnidt(2:Zh+1) = dnidt(2:Zh+1)+DCX_ip2_i-DCX_i_im2;
            %Taking double electron charge exchange only for HCI
            
            Zl = obj.HCI.Zl;
            id = Zh+1;
            EII_im1_i_l = obj.EII_rate((1:Zl)+id);
            EII_i_ip1_l = obj.EII_rate((2:(Zl+1))+id);
            RR_ip1_i_l = [obj.RR_rate((3:(Zl+1))+id),0];
            RR_i_im1_l = obj.RR_rate((2:(Zl+1))+id);
            CX_ip1_i_l = [obj.CX_rate((3:(Zl+1))+id),0];
            CX_i_im1_l = obj.CX_rate((2:(Zl+1))+id);
            ESC_a_l = obj.ESC_a_rate((2:(Zl+1))+id);
            ESC_r_l = obj.ESC_r_rate((2:(Zl+1))+id);
            % dnidt(Zh+2) = -obj.EII_rate(1+id)+obj.RR_rate(2+id)+obj.CX_rate(2+id);
            dnidt(Zh+2) = 0;

            dnidt(Zh+3:end) = EII_im1_i_l-EII_i_ip1_l+RR_ip1_i_l-RR_i_im1_l+CX_ip1_i_l-CX_i_im1_l+ESC_a_l+ESC_r_l;
            DCX_ip2_i = [obj.DCX_rate((4:Zl+1)+id),0,0];
            DCX_i_im2 = [0,obj.DCX_rate((3:Zl+1)+id)];

            dnidt(Zh+3:end) = dnidt(Zh+3:end)+DCX_ip2_i-DCX_i_im2;
        end

    end
end