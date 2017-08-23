%This function lists the A,B,C,n coefficients
%required to solve for the density of a organic liquid.
%The liquid is: Hexadecane (C16H34)

%NOTE: The temperature range of validity of these coefficients/formula
%      is 291.31 - 723.00 K

function [A,B,C,n,Tlow,Thigh]= HexRhoCoefs()

%Relevant coefficients to calculate density

Tlow = 291.31;
Thigh = 723.00;
A = 0.2471;
B = 0.26626;
C = 723.00;
n = 0.28571;

