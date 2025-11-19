classdef ElectronBeamClass < handle
    properties
        %Parameter:
        %Energy (eV) electron energy
        %rH (m) electron beam radius using Herrman Radius formula
        %Velocity (m/s) the velocity of the electron
        %B_filed (T) The magentic field in the trap center
        %I_current (A) the current of the electron beam
        %Tc (K) the temperature of the cathode
        %rc (m) the Radius of the cathode emision surface
        %Bc (T) the magnetic field at the cathode position
        %je (A/cm^2) the current density in the central trap paet
        %jc (A/cm^2) the current density in the cathode part
        %ne (/cm^3) the electron density 
        Energy = 200e+3;
        rH;
        Velocity;
        B_field = 6;
        I_current=200e-3;
        Tc=1350;
        rc= 1.7e-3;%Baumann thesis
        Bc=1e-3;
        je;
        jc;
        ne;
        me_kg = 9.109e-31;
        me_eV ;
        c= 299792458;
        e=1.6022e-19;
        epsilon0=8.854e-12;
        kb = 1.381e-23;
    end
    methods
        function obj= InitialParam(obj)
            %Convert electron mass from kg to eV
            obj.me_eV=obj.me_kg*obj.c.^2./obj.e;
            %Calculate the electron velocity using rel version formula
            obj.Velocity=rel_v();
            rB=BrillouinRadius();
            obj.rH=HerrmanRadius(rB);
            obj.ne=ElectronDensity();
            obj.jc=ElectronDensityCathode();

            function velocity=rel_v()
                %Calculate the electron velocity using rel version formula
                E_total=obj.Energy+obj.me_eV;
                beta=sqrt(1-(obj.me_eV/E_total)^2);
                velocity=beta*obj.c;
            end

            function rB=BrillouinRadius()
                %The estimation of the electron beam radius using Brillouin
                %Radius fromula
                rB=sqrt(2*obj.me_kg*obj.I_current/pi...
                    /obj.epsilon0/obj.e/obj.Velocity/obj.B_field.^2);
            end
            function rH=HerrmanRadius(rB)
                %The estimation of the electron beam radius using Herrman
                %Radius fromula
                %Ref: https://doi.org/10.1063/1.5026961
                rH=rB*sqrt(0.5+sqrt(0.25+8*obj.me_kg*obj.kb*obj.Tc*obj.rc^2/(obj.e*obj.B_field* ...
                    rB^2)^2+obj.Bc^2*obj.rc^4/obj.B_field^2/rB^4));
            end

            function ne=ElectronDensity()
                %Calculation for the current area density and electron density
                obj.je=obj.I_current./(pi*obj.rH.^2);
                ne=obj.je/obj.e/obj.Velocity;
                obj.je=obj.je/1e+4; %convert from A/m^2 to A/cm^2
                ne=ne/1e+6;%Convert from /m^3 to /cm^3
            end

            function jc=ElectronDensityCathode()
                jc=obj.I_current./(pi*obj.rc.^2);
                jc=jc/1e+4; %convert from A/m^2 to A/cm^2
            end
        end
    end
end