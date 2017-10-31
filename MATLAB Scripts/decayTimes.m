clear
clc
format compact
format long

To = 298;       %initial chamber temperature [K]

%initial droplet diameter [mm]
do = [2.215
1.752
2.096
3.071
3.181
2.969
2.762
2.091
4.136
3.821
2.261
2.411
4.179
4.295
1.908
2.182
1.727
3.731
2.116];

%burning rate before flame contraction [mm^2/s]
K = [0.857907745
0.930168043
0.66226582
0.984639004
1.198449199
1.06030366
0.503190516
0.614689639
0.941413661
0.648033754
0.622164639
0.815213201
1.064993191
0.995278088
0.335056535
0.476306483
0.329643425
0.664302482
0.284856146];

%specify mass fraction of solute (Hexadecane) in mixture
sizeBlock1 = 8;
for i=1:length(do)
    if i <= sizeBlock1
        YsuVector(i) = 0.05;      
    else
        YsuVector(i) = 0.20;
    end
end


%import in molecular weights of each compoents and calculate mole fractions
%molecular weights in units of kg/kmol
[PcSU, TcSU,MWsu] = HexCriticalProps();
[PcSV, TcSV,MWsv] = HepCriticalProps();

%import in coefficients to calculate densities of each component
%at reference temperature T=T0=298K
T0 = 298;
[ArhoSV, BrhoSV, CrhoSV, nrhoSV, TrhoLSV, TrhoHSV ] = HepRhoCoefs();
%compute Solvent (Heptane) density at 298 K [kg/m^3}
rhoSvT0 = rhoCalc(ArhoSV, BrhoSV, CrhoSV,...
    nrhoSV, TrhoLSV, TrhoHSV, T0)*1000;
[Asu,Bsu,Csu,nsu,TrhosuL,TrhosuH]= HexRhoCoefs() ;
%compute Solute (Hexadecane) density at 298 K [kg/m^3]
rhoSuT0 = rhoCalc(Asu,Bsu,Csu,...
    nsu,TrhosuL,TrhosuH, T0)*1000;


%import in Andrade coefficients and calculate liquid viscosity
%if each component in (kg/m-s) at ref temp T=T0
%(1) => solute (hexadecane) and (2) => solvent (heptane)
[AandSU,BandSU,TminandSU,TmaxandSU]= HexAndradeCoefs();
[mu1o ]= AndradeVSL(AandSU,BandSU,TminandSU,TmaxandSU, T0)*.001; %N-s/m^2
[AandSV,BandSV,TminandSV,TmaxandSV]= HepAndradeCoefs();
[mu2o ]= AndradeVSL(AandSV,BandSV,TminandSV,TmaxandSV, T0)*.001; %N-s/m^2

for i=1:length(YsuVector)
    
    Ysu = YsuVector(i);
    
    %%%%%%%%%%%%%%%%%% calculate mole fractions %%%%%%%%%%%%%%%%%%%%%
    MWmix = 1/( (Ysu/MWsu) + ((1-Ysu)/MWsv) );
    X1 = (Ysu*MWmix/MWsu);
    
    %%%%%%%%%%%%%%%%%% calculate density of mixture %%%%%%%%%%%%%%%%%%%%%
    
    %mixture density [kg/m^3]
    rhoMixT0 = 1/ ( ((1-Ysu)/rhoSvT0) + (Ysu/rhoSuT0) )
    
    %%%%%%%%%%%%%%%% calculate mixture dynamic viscosity %%%%%%%%%%%%%%%%
    
    %calculate mixture static viscosity in kg/m-s or N-s/m^2
    muMixture = exp( X1*log(mu1o) + (1-X1)*log(mu2o) + ...
        X1*(1-X1)*1.08*(1.343-X1*0.685) );
    
    %compute kinematic viscosity in mm^2/s
    nuMix(i) = ( muMixture/rhoMixT0 )*(1e6);
    
end

Uo = 800;    %initial characteristic velocity in droplet [mm/s]

%Peclet number := Udo/K =convection/diffusive effects
%Peclet no small => Convection effects are small. We want
%beta to be small.  This is a number we impose.
beta = 0.1;  

for i = 1:length(do)
    tv(i) = ( do(i)^2/(4*nuMix(i)) )*log( Uo*do(i) / ( 2*K(i)*beta ) );
end

disp('viscous decay times are (seconds): ')
disp(tv')



































