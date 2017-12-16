%This function computes the saturated temperature, Tsat, given 
%the saturated pressure as input using the Antoinne relation (KDB)

%ln(Psat) = A*ln(T) + B/T + C + D*T^2 where Pvp in kPa, T in K.

%The KDB correlation calculates the saturated pressure, Pvp given
%the temperature, Tsat.  However, in this case, we wish to calculate
%T=Tsat, by fixing Pvp= pressure of the chamber.  The equation is
%non-linear in T and is solved using Newton's Method with damping.

%INPUT: Psat - Saturated Pressure [kPa]
%       A,B,C,D, Tlow, Thigh - Correlation coefficients in Antoinne Equation       

function [Tsat] = AntoinneTsat(Psat,A,B,C,D,Tlow,Thigh)

%using Newton's Method to determine saturated temperature
f = @(T) A*log(T) + B/T + C + D*T^2 - log(Psat);
Df = @(T) A/T - B/(T*T) + 2*D*T;

%tolerance parameters input into Newton's Method
FTOL    = eps;     %here we'll tak FTOL as the max norm
XTOL    = eps;     %here we use the max norm
kmax    = 100;     %max iteration
lmax    = 15;      %minimum damping factor
T0 = (Thigh+Tlow)/2;    %initial guess for Tsat as input into Newton's Method

[Tsat, k, EuNormFx, DX, dampfactor] = ...
    dampNewton(T0,Df,f,FTOL,XTOL,kmax,lmax);
if (Tsat < Thigh) && (Tsat > Tlow)
    %do nothing
else
    disp('Solvent Tsat is out of range of validity of the Antoinne Equation')
end