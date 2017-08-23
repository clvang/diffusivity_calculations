%This function lists the A,B,C,n coefficients
%required to solve for the density of a organic liquid.
%The liquid is: N-Heptane (C7H16)

%NOTE: The temperature range of validity of these coefficients/formula
%      is [153.75-535.20 K]

function [A,B,C,n,Tlow,Thigh]= HepRhoCoefs()

%Relevant coefficients to calculate density
Tlow = 153.75;
Thigh = 535.20;
A = 0.2341;
B = 0.26058;
C = 540.2;
n = 0.28571;



