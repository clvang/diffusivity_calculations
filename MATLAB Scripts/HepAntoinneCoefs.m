%This function lists the A,B,C,D coefficients
%required to solve for the saturated Pressure or,
%saturated Temperature of an organic liquid, using Antoinne Equation
%The liquid is: N-Heptane (C7H16)

%NOTE: The temperature range of validity of these coefficients/formula
%      is 182.56-540.26 K

function [A,B,C,D, Tlow, Thigh]= HepAntoinneCoefs()

%TEMP RANGE OF VALIDITY: 182.56-540.26 K
Tlow = 182.56;
Thigh = 540.26;

A = -1.412388E+01;
B = -8.030070E+03;
C = 1.081461E+02;
D = 1.204855E-05;