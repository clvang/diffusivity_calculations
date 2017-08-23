%This funciton calculates the liquid density of an Organic liquid

%INPUT:  
%   T - Temperature at which to the density is to be evaluated [Kelvin]
%   A,B,C,n,Tlow,Thigh - correlation coefficients for density computation

%OUTPUT: density, rho [g/ml] 
%NOTE:  g/ml = 1000 kg/m^3

function     [rho] = rhoCalc(A,B,C,n,Tlow,Thigh,T)

if (T < Thigh) && (T > Tlow)
    rho = A*(B^(-(1-T/C)^n));       %g/ml
else
    disp('T is out of temp. range of validity for Heptane')
end
