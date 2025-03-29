% analytic solution for FvCB model 
% first drafted by ZHT 
% Last update by ZHT Mar 26, 2025
% Copyright reserved by ZHT @ YNU
function [A1,A2,A3] = fAnalytic_BBerry_Vc(g0, g1, gb, gm, Rd, C, V, G, rh, D, Ca);
% INPUT ******************************************************************
% g0: intercept of stomatal model                                        *
% g1: slope of stomatal model                                            *
% gb: boundary layer conductance (mol m-2 s-1)                           *
% gm: mesophyll conductance (mol m-2 s-1)                                *
% Rd: respiration (umol m-2 s-1)                                         *
% C: Cv (ubar) very sensitive parameter
% V: Vcmax (umol m-2 s-1)
% G: gamma_str (ubar)
% rh: relative humidity (range 0 and 1) drive Ball-Berry stomatal model
% D: water vapor pressure deficit (kPa)
% Ca: ambient CO2 concentration (ubar)
% OUTPUT *****************************************************************
% A1-A3: Three solutions for the Qubic equation                          *
kk = (g0*gb + g0*gm + gb*gm - g1*gb^2*rh - g1*gb*gm*rh); 
p = -(V*g0*gb - Rd*g0*gm - Rd*gb*gm - Rd*g0*gb + V*g0*gm + V*gb*gm + Ca*g0*gb^2 + Ca*gb^2*gm + Rd*g1*gb^2*rh - V*g1*gb^2*rh + C*g0*gb*gm + 2*Ca*g0*gb*gm + Rd*g1*gb*gm*rh - V*g1*gb*gm*rh - C*g1*gb^2*gm*rh - Ca*g1*gb^2*gm*rh)/kk; 
q = (gb*(Ca^2*g0*gb*gm - C*Rd*g0*gm - Ca*Rd*g0*gb - 2*Ca*Rd*g0*gm - Ca*Rd*gb*gm + Ca*V*g0*gb + 2*Ca*V*g0*gm + Ca*V*gb*gm - G*V*g0*gm + C*Ca*g0*gb*gm + C*Rd*g1*gb*gm*rh + Ca*Rd*g1*gb*gm*rh - Ca*V*g1*gb*gm*rh + G*V*g1*gb*gm*rh))/kk; 
r = (Ca*g0*gb^2*gm*(C*Rd + Ca*Rd - Ca*V + G*V))/kk; 

% 计算判别式和根
Q = (p^2 - 3*q) / 9; 
R = (2*p^3 - 9*p*q + 27*r)/54;
TH = acos(R/sqrt(Q^3)); 

A1 = -2*sqrt(Q)*cos(TH/3) - p/3; 
A2 = -2*sqrt(Q)*cos((TH + 2*pi)/3) - p/3; 
A3 = -2*sqrt(Q)*cos((TH + 4*pi)/3) - p/3; 
