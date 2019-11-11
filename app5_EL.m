close all
clear all
clc
%% Script de l'app 5 S5e
% Par Hubert Dube
% Debute le 7/11/2019
specs_app5;
%% Telescope A ------------------------------------------------------------
%--------------------------------------------------------------------------
%-----------------------------ELEVATION------------------------------------
%--------------------------------------------------------------------------
%% traduction des specifications
trad_specs;

% % creation du lieu des racines
% figure(1)
% hold on
% rlocus(FTBO_EL)
% scatter(real(s),imag(s),'p')
% title('LdR FTBO-EL')
%% Compensateur pour ELEVATION
% creation d'un avance de phase
marge = 15;%20.5;%19;5
ka_red = 0;
phase_EL = rad2deg(angle(numEL/polyval(denEL,s(1))));
delta_phi_AvPh_EL = -180 - phase_EL + 360 + marge;
phi_AvPh_EL = 180 - rad2deg(atan2(imag(s(1)),real(s(1))));
alpha_AvPh_EL = 180-phi_AvPh_EL;
phi_z_AvPh_EL = (alpha_AvPh_EL + delta_phi_AvPh_EL)/2;
phi_p_AvPh_EL = (alpha_AvPh_EL - delta_phi_AvPh_EL)/2;
z_AvPh_EL = real(s(1)) - imag(s(1))/tan(deg2rad(phi_z_AvPh_EL));
p_AvPh_EL = real(s(1)) - imag(s(1))/tan(deg2rad(phi_p_AvPh_EL));


ka_AvPh_EL = 1/norm((s(1)-z_AvPh_EL)/(s(1)-p_AvPh_EL)* numAZ/polyval(denAZ,s(1)));
num_AvPh_EL = [1 -z_AvPh_EL];
den_AvPh_EL = [1 -p_AvPh_EL];
AvPh_EL = ka_AvPh_EL*tf(num_AvPh_EL,den_AvPh_EL);

% % creation du lieu des racines avec l'AvPH
% figure(2)
% hold on
% rlocus(FTBO_EL*AvPh_EL)
% scatter(real(s),imag(s),'p')
% title('LdR FTBO-EL*AvPh-EL')

% figure(3)
% hold on
% step(feedback(FTBO_EL*AvPh_EL,1))
% step(feedback(FTBO_EL,1))
% legend('FTBO-EL*AvPh','FTBO-EL')
% title('Step comparaison FTBO-EL*AvPh-EL')
%% verification du systeme asservi par l'AvPh
FTBF_EL_AvPh = feedback(FTBO_EL*AvPh_EL,1);
stepinfo(FTBF_EL_AvPh)

[num_FTBO_AvPh_EL,den_FTBO_AvPh_EL] = tfdata(FTBO_EL*AvPh_EL,'v');
kvel_EL = polyval(num_FTBO_AvPh_EL,0)/polyval([den_FTBO_AvPh_EL(1:end-1)],0);
eru_EL = 1/kvel_EL;
% on conclu qu l'ordre de la TF doit etre augmente pcq il y a une erreur en
% reponse  a une rampe unitaire

%% creation d'un PI
z_PI_EL = real(s(1))/(6)%(6.35)
ka_PI_EL = 1/norm((s(1)-z_PI_EL)/(s(1))* polyval(num_FTBO_AvPh_EL,s(1))/polyval([den_FTBO_AvPh_EL],s(1)))
PI_EL = ka_PI_EL* tf([1 -z_PI_EL],[1 0])

%% ajout d'un coupe bande
% trouver la frequence de coupure
freq_coup = 123 % rad/sec trouve avec bode
w_width = 40
num_band_stop = [1 0 freq_coup^2];
den_band_stop = [1 w_width freq_coup^2];
band_stop = tf(num_band_stop,den_band_stop);

figure(4);
hold on;
step(feedback(FTBO_EL*AvPh_EL*PI_EL,1),5);
step(feedback(FTBO_EL*AvPh_EL*PI_EL*band_stop,1),5);
legend('PI_EL','band_stop');
stepinfo(feedback(FTBO_EL*AvPh_EL*PI_EL*band_stop,1))

%% verification du systeme asservi par l'AvPh et un PI
FTBF_EL_PI = feedback(FTBO_EL*AvPh_EL*PI_EL*band_stop,1)
stepinfo(FTBF_EL_PI)

[num_FTBO_PI_EL,den_FTBO_PI_EL] = tfdata(FTBO_EL*AvPh_EL*PI_EL*band_stop,'v');
kacc = polyval(num_FTBO_PI_EL,0)/polyval([den_FTBO_PI_EL(1:end-2)],0);
epu = 1/kacc

figure(6)
margin(FTBO_EL*AvPh_EL*PI_EL*band_stop)
[Gm2_EL,Pm2_EL,Wp2_EL,Wg2_EL] = margin(FTBO_EL*AvPh_EL*PI_EL*band_stop);
RM_EL = Pm2_EL/Wg2_EL*pi/180


ramp = [0:0.1:20];
figure
hold on x
lsim(FTBF_EL_PI,ramp,ramp)


y_ramp = lsim(FTBF_EL_PI,ramp,ramp);
y_ramp_diff = y_ramp-ramp';
para = 0.5 * ramp'.^2;
y_para = lsim(FTBF_EL_PI,para,ramp);
y_para_diff = y_para-para;

figure()
hold on
plot(ramp',y_ramp_diff)
plot(ramp',y_para_diff)