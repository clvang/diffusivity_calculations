%This function lists the A,B,C,D coefficients
%required to solve for the saturated Pressure or,
%saturated Temperature of an organic liquid, using Antoinne Equation
%The liquid is: Hexadecane (C16H34)

%NOTE: The temperature range of validity of these coefficients/formula
%      is 291.32 K - 720.60 K

function [A,B,C,D, Tlow, Thigh]= HexAntoinneCoefs()


%TEMP RANGE OF VALIDITY: 291.32 K - 720.60 K
Tlow = 291.32;
Thigh = 720.60;

A = -2.059165E+01;
B = -1.552894E+04;
C = 1.601365E+02;
D = 8.014932E-06;