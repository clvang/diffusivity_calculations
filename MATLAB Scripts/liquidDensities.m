clear
clc
format compact
format long

%computes ratio of (rho_0/rho_sat)^(1/3) to account for increases in do due
%to heating

P   = [1, 3];     %chamber pressures [atm]
T0 = 298;        %initial chamber temperature [K]

% % for heptane 80, hexadecane 20
% Ysv = 0.80;     %solvent mass fraction (Heptane)
% Ysu = 0.20;      %solute mass fraction (Hexadecane)

% for heptane 95, hexadecane 05
Ysv = 0.80;     %solvent mass fraction (Heptane)
Ysu = 0.20;      %solute mass fraction (Hexadecane)

for i = 1:length(P)
    
    disp('*******************************************************')
    fprintf('The pressure is: %d atm \n',P(i))

    %% SOLVENT PROPERTIES
    
    
    %Properties of more volatile (usually higher mass fraction or SOLVENT)
    %component.
    %COMPONENT NAME: Heptane (C7H16)
     
    %import in A,B,C,n,Tlow,Thigh correlation coefficients for density
    [ArhoSV, BrhoSV, CrhoSV, nrhoSV, TrhoLSV, TrhoHSV ] = HepRhoCoefs();
    
    rhoSvT0 = rhoCalc(ArhoSV, BrhoSV, CrhoSV,...  %g/ml
            nrhoSV, TrhoLSV, TrhoHSV, T0);   %compute Solvent density at T0
    
    %import in correlation coefficients in Anttoine coeffcients
    [AsatSV, BsatSV, CsatSV, DsatSV, TpsatsvL, TpsatsvH]=...
        HepAntoinneCoefs();
    
    %set Psat=chamber pressure
    Pvp = P(i)*101.325;   %specify saturated pressure, [kPa]
    
    %compute Tsat using Antoinne equation
    TsatSV = AntoinneTsat(Pvp,AsatSV,BsatSV,CsatSV,...
            DsatSV, TpsatsvL, TpsatsvH);
   
    %compute rho_sat using Tsat obtained above
    rhoSvTsat = rhoCalc(ArhoSV, BrhoSV, CrhoSV,...
            nrhoSV, TrhoLSV, TrhoHSV, TsatSV);      %g/ml

    %% SOLUTE PROPERTIES
    
    %Properties of less volatile (usually lower mass fraction or SOLUTE)
    %component.
    %COMPONENT NAME: Hexadecane (C16H34)
       
    [Asu,Bsu,Csu,nsu,TrhosuL,TrhosuH]= HexRhoCoefs() ;
    
    rhoSuT0 = rhoCalc(Asu,Bsu,Csu,...   %g/ml
            nsu,TrhosuL,TrhosuH, T0);   %compute Solute density at T0
        
    [AsatSU,BsatSU,CsatSU,DsatSU, TpsatsuL, TpsatsuH]= HexAntoinneCoefs();
    
    TsatSU = AntoinneTsat(Pvp,AsatSU,BsatSU,CsatSU,...
            DsatSU, TpsatsuL, TpsatsuH);
    
        
    rhoSuTsat = rhoCalc(Asu,Bsu,Csu,...     %g/ml
            nsu,TrhosuL,TrhosuH, TsatSU);   %compute Solute density at Tsat

    
    %% MIXTURE PROPERTIES
    
    rhoMixT0 = 1/ ( (Ysv/rhoSvT0) + (Ysu/rhoSuT0) );        %g/ml
    rhoMixTsat = 1/( (Ysv/rhoSvTsat) + (Ysu/rhoSuTsat) );   %g/ml
    
    
    %% Compute percent increase in do due to droplet heating
    
    rhoOrhoSat = (rhoSvT0 / rhoSvTsat)^(1/3);
    percentRhoOrhoSat = (rhoOrhoSat-1)*100;
    disp('% d0 increase, considering just the Solvent is: ')
    disp(percentRhoOrhoSat)
    
    rhoOrhoSatMixture = (rhoMixT0/rhoMixTsat)^(1/3);
    percentrhoOrhoSatMixture = (rhoOrhoSatMixture-1)*100;
    disp('% d0 increase, considering the MIXTURE is: ')
    disp(percentrhoOrhoSatMixture)
    disp('*******************************************************')
    
end







