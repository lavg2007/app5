close all
clear all
clc
%% Script de l'app 5 S5e
% Par Hubert Dube
% Debute le 7/11/2019
specs_app5
%% Telescope A ------------------------------------------------------------
%--------------------------------------------------------------------------
%-----------------------------ELEVATION------------------------------------
%--------------------------------------------------------------------------
%% traduction des specifications
trad_specs

% creation du lieu des racines
figure(1)
hold on
rlocus(FTBO_EL)
scatter(real(s),imag(s),'p')
title('LdR FTBO-EL')
%% Compensateur pour ELEVATION
% creation d'un avance de phase
phase_EL = rad2deg(angle(numEL/polyval(denEL,s(1))))
delta_phi_AvPh_EL = -180 - phase_EL +360
phi_AvPh_EL = 180 - rad2deg(atan2(imag(s(1)),real(s(1))));
alpha_AvPh_EL = 180-phi_AvPh_EL;
phi_z_AvPh_EL = (alpha_AvPh_EL + delta_phi_AvPh_EL)/2;
phi_p_AvPh_EL = (alpha_AvPh_EL - delta_phi_AvPh_EL)/2;
z_AvPh_EL = real(s(1)) - imag(s(1))/tan(deg2rad(phi_z_AvPh_EL));
p_AvPh_EL = real(s(1)) - imag(s(1))/tan(deg2rad(phi_p_AvPh_EL));


ka_AvPh_EL = 1/norm((s(1)-z_AvPh_EL)/(s(1)-p_AvPh_EL)* numAZ/polyval(denAZ,s(1)))
num_AvPh_EL = [1 -z_AvPh_EL];
den_AvPh_EL = [1 -p_AvPh_EL];
AvPh_EL = ka_AvPh_EL*tf(num_AvPh_EL,den_AvPh_EL);

% creation du lieu des racines avec l'AvPH
figure(2)
hold on
rlocus(FTBO_EL*AvPh_EL)
scatter(real(s),imag(s),'p')
title('LdR FTBO-EL*AvPh-EL')

figure(3)
hold on
step(feedback(FTBO_EL*AvPh_EL,1))
step(feedback(FTBO_EL,1))
legend('FTBO-EL*AvPh','FTBO-EL')
title('Step comparaison FTBO-EL*AvPh-EL')
%% verification du systeme asservi par l'AvPh
FTBF_EL_AvPh = feedback(FTBO_EL*AvPh_EL,1)
stepinfo(FTBF_EL_AvPh)
