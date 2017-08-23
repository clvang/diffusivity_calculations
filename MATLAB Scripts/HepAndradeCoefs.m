%This function lists the thermodynamic properties of Heptane
%required to calculate the viscosity of Heptane using the Andrade
%Equation



function [A,B,Tmin,Tmax]= HepAndradeCoefs()

%Coefficients and temp range from CHERIC
Tmin = 183.15;      %minmum temp of validity, K
Tmax = 373.15;      %max temp of validity, K

A = -4.325;
B = 1006;


