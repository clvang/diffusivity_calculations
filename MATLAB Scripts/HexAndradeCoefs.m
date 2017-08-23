%This function lists the thermodynamic properties of Hexadecane
%required to calculate the viscosity of Heptane using the Andrade
%Equation



function [A,B,Tmin,Tmax]= HexAndradeCoefs()

%Coefficients and temp range from CHERIC
Tmin = 293.15;      %minmum temp of validity, K
Tmax = 553.15;      %max temp of validity, K

A = -4.643;
B = 1700;


