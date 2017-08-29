clear
clc
format compact
format long


%this script calculates the liquid diffusion coefficients
%at infinite dilution.

P   = [1];     %chamber pressures [atm]


for i = 1:length(P)
    %% SOLVENT PROPERTIES (B) (Heptane)
    
    
    %Properties of more volatile (usually higher mass fraction or SOLVENT)
    %component.
    %COMPONENT NAME: Heptane (C7H16)
    
    %import in A,B,C,n,Tlow,Thigh correlation coefficients for density
    %These coefs are needed in "rhoCalc" function below
    [ArhoSV, BrhoSV, CrhoSV, nrhoSV, TrhoLSV,...
        TrhoHSV] = HepRhoCoefs();
    
    %import in correlation coefficients in Anttoine coeffcients.
    %compute Tsat using Antoinne equation and compute corresponding
    %density in kg/m^3 at Tsat. Tsat is the NORMAL BP
    %e.g. when the pressure is 1 atm=101.325 kPa
    [AsatSV, BsatSV, CsatSV, DsatSV, TpsatsvL, TpsatsvH]=...
        HepAntoinneCoefs();
    TsatSV = AntoinneTsat(101.325,AsatSV,BsatSV,CsatSV,...
        DsatSV, TpsatsvL, TpsatsvH);    %Normal BP temp [K] at 101.325 kPa
    %NOTE: Vargaftik, TsatSV = 98.43 deg C = 371.58 K (CHECK)
    
    
    rhoSvTsat = rhoCalc(ArhoSV, BrhoSV, CrhoSV,...
        nrhoSV, TrhoLSV, TrhoHSV, TsatSV)*1000; %density [kg/m^3] at normal BP
    %NOTE: Vargaftik, rhoSvTsat = 6.139e2 kg/m^3 (CHECK)
    
    %Compute boiling point Tbp corresponding to the champer pressure
    %using the Antoinne equation
    TbpSVChamber = AntoinneTsat(P(i)*101.325,AsatSV,BsatSV,CsatSV,...
        DsatSV, TpsatsvL, TpsatsvH);    %BP temp [K] at chamber pressure

    
    %import in Andrade coefficients and calculate liquid viscosity
    %in centiPoise (cp) of the SOLVENT at T=TbpSVchamber NOTE: Here we assume
    %that T=TbpSVchamber.  When the droplet is already heated, T should be at the
    %boiling point temperature corresponding to the chamber pressure
    [AandSV,BandSV,TminandSV,TmaxandSV]= HepAndradeCoefs();
    [VSLsv ]= AndradeVSL(AandSV,BandSV,TminandSV,TmaxandSV, TbpSVChamber); %[cP]
    %NOTE: Vargaftik, VSLsv = 0.201 cP (CHECK)
    
    
    %calculate surface tension using corresponding states of Bird and Brock
    [PcSV, TcSV,MWsv] = HepCriticalProps();
    sigmaSV = 12.6191; %BrockBirdCS(PcSV,TcSV,TbpSVChamber);  %[erg/cm^2] @ normal BP
    %NOTE: Vargaftik, sigmaSV = 12.6191e-3 N/m = 12.6191 erg/cm^2 (CHECK)
    
    
    
    %% SOLUTE PROPERTIES  (A) (Hexadecane)
    
    %Properties of less volatile (usually lower mass fraction or SOLUTE)
    %component.
    %COMPONENT NAME: Hexadecane (C16H34)
    
    [Asu,Bsu,Csu,nsu,TrhosuL,TrhosuH]= HexRhoCoefs() ;
    
    [AsatSU,BsatSU,CsatSU,DsatSU, TpsatsuL, TpsatsuH]= HexAntoinneCoefs();
    
%     TsatSU = AntoinneTsat(101.325,AsatSU,BsatSU,CsatSU,...
%         DsatSU, TpsatsuL, TpsatsuH);        
    %NOTE: Vargaftik, TsatSU = 287.05 deg C = 560.2 K (CHECK)
    
    
    rhoSuTsat = rhoCalc(Asu,Bsu,Csu,...     %kg/m^3
        nsu,TrhosuL,TrhosuH, TbpSVChamber)*1000;   %compute Solute density at Tsat
    %NOTE: Vargaftik, rhoSuTsat = 7.19e2 kg/m^3 (CHECK)
    
%     %Compute boiling point Tbp corresponding to the champer pressure
%     %using the Antoinne equation
%     TbpSUChamber = AntoinneTsat(P(i)*101.325,AsatSU,BsatSU,CsatSU,...
%         DsatSU, TpsatsuL, TpsatsuH);    %BP temp [K] at chamber pressure
    
    %import in Andrade coefficients and calculate liquid viscosity
    %in centiPoise (cp) of the SOLUTE at T=TbpSUchamber NOTE: Here we assume
    %that T=TbpSVchamber.  When the droplet is already heated, T should be at the
    %boiling point temperature corresponding to the chamber pressure.
    %We evaluate VSLsu at the boiling point of the Solvent (Heptane)
    [AandSU,BandSU,TminandSU,TmaxandSU]= HexAndradeCoefs();
    [VSLsu ]= AndradeVSL(AandSU,BandSU,TminandSU,TmaxandSU, TbpSVChamber); %[cP]
    %NOTE: Vargaftik, VSLsu = 0.911e-3 N-s/m^2 = 0.911 cP (CHECK)
    
    %calculate surface tension using corresponding states of Bird and Brock
    [PcSU, TcSU,MWsu] = HexCriticalProps();

    sigmaSU = 21.03; %BrockBirdCS(PcSU,TcSU,TbpSVChamber)  %[erg/cm^2] @ normal BP
    %NOTE: Vargaftik, sigmaSU=21.03e-3 N/m = 21.03 erg/cm^2
    
    
% DAB = tynCalus(MWsolv, rhoSATsolv,MWsolu,rhoSATsolu,...
%     sigmaSolv,sigmaSolu,VSLsv,T)

    %find D of Solute A (hexadecane) diffusing into Solvent B (Heptane) in
    %cm^2/s
    DAB = tynCalus(MWsv, rhoSvTsat,MWsu,rhoSuTsat,...
        sigmaSV,sigmaSU,VSLsv,TbpSVChamber);
    DAB = DAB/(10^4); %m^2/s
    disp('The molecular diffusivity D_{AB} [m^2/s] is: ')
    disp(DAB)
    
    %find D of Solute B (Heptane) diffusing into Solvent A (Hexadecane) in
    %cm^2/s
    DBA = tynCalus(MWsu, rhoSuTsat,MWsv,rhoSvTsat,...
        sigmaSU,sigmaSV,VSLsu,TbpSVChamber);
    DBA = DBA/(10^4);
    disp('The molecular diffusivity D_{BA} [m^2/s} is: ')
    disp(DBA)
    
end



% Vb = (MWsolv/rhoSATsolv)/(1e-3); %molar volume of SOLVENT cm^3/mol
%                                  %at normal BP
% 
% Va = (MWsolu/rhoSATsolu)/(1e-3); %molar volume of SOLUTE cm^3/mol
%                                  %at normal BP
% 
