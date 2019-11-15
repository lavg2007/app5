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
Profile_Tracking;

% % creation du lieu des racines
% figure(1)
% hold on
% rlocus(FTBO_EL)
% scatter(real(s),imag(s),'p')
% title('LdR FTBO-EL')
figure()
step(feedback(FTBO_EL,1))
saveas(gcf,'FTBO_EL.png')
figure()
margin(FTBO_EL)
saveas(gcf,'marges_FTBO_EL.png')

%% Compensateur pour ELEVATION
% creation d'un avance de phase
marge = 20%20.5;%19;5
phase_EL = rad2deg(angle(numEL/polyval(denEL,s(1))));
delta_phi_AvPh_EL = -180 - phase_EL + 360 + marge;
phi_AvPh_EL = 180 - rad2deg(atan2(imag(s(1)),real(s(1))));
alpha_AvPh_EL = 180-phi_AvPh_EL;
phi_z_AvPh_EL = (alpha_AvPh_EL + delta_phi_AvPh_EL)/2;
phi_p_AvPh_EL = (alpha_AvPh_EL - delta_phi_AvPh_EL)/2;
z_AvPh_EL = real(s(1)) - imag(s(1))/tan(deg2rad(phi_z_AvPh_EL));
p_AvPh_EL = real(s(1)) - imag(s(1))/tan(deg2rad(phi_p_AvPh_EL));


ka_AvPh_EL = 1/norm((s(1)-z_AvPh_EL)/(s(1)-p_AvPh_EL)* numEL/polyval(denEL,s(1)));
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

figure(30)
hold on
margin(FTBO_EL*AvPh_EL)
saveas(gcf,'FTBO_EL_AvPh_marge_no.png')

[num_FTBO_AvPh_EL,den_FTBO_AvPh_EL] = tfdata(FTBO_EL*AvPh_EL,'v');
kvel_EL = polyval(num_FTBO_AvPh_EL,0)/polyval([den_FTBO_AvPh_EL(1:end-1)],0);
eru_EL = 1/kvel_EL;
% on conclu qu l'ordre de la TF doit etre augmente pcq il y a une erreur en
% reponse  a une rampe unitaire

%% creation d'un PI
z_PI_EL = real(s(1))/(4.4)%(6.35)
ka_PI_EL = 1/norm((s(1)-z_PI_EL)/(s(1))* polyval(num_FTBO_AvPh_EL,s(1))/polyval([den_FTBO_AvPh_EL],s(1)))
PI_EL = ka_PI_EL* tf([1 -z_PI_EL],[1 0])

%% ajout d'un coupe bande
% trouver la frequence de coupure
freq_coup = 123 % rad/sec trouve avec bode
w_width = 40
num_band_stop = [1 0 freq_coup^2];
den_band_stop = [1 w_width freq_coup^2];
band_stop = tf(num_band_stop,den_band_stop);

% [num,den]= cheby1(2,0,[freq_coup-w_width/2 freq_coup+w_width/2],'stop','s');
%  band_stop = tf(num,den);

figure(4);
hold on;
step(feedback(FTBO_EL,1),5);
step(feedback(FTBO_EL*AvPh_EL,1),5);
step(feedback(FTBO_EL*AvPh_EL*PI_EL,1),5);
step(feedback(FTBO_EL*AvPh_EL*PI_EL*band_stop,1),5);
legend('FTBO_EL','AvPh_EL','PI_EL','band_stop');
xlim([0 3])
saveas(gcf,'step_A_EL_comp.png')
stepinfo(feedback(FTBO_EL*AvPh_EL*PI_EL*band_stop,1))
%% verification du systeme asservi par l'AvPh et un PI
FTBF_EL_PI = feedback(FTBO_EL*AvPh_EL*PI_EL*band_stop,1)
stepinfo(FTBF_EL_PI)

[num_FTBO_PI_EL,den_FTBO_PI_EL] = tfdata(FTBO_EL*AvPh_EL*PI_EL*band_stop,'v');
kacc = polyval(num_FTBO_PI_EL,0)/polyval([den_FTBO_PI_EL(1:end-2)],0);
epu = 1/kacc
% 
figure(6)
margin(FTBO_EL*AvPh_EL*PI_EL*band_stop)
[Gm2_EL,Pm2_EL,Wp2_EL,Wg2_EL] = margin(FTBO_EL*AvPh_EL*PI_EL*band_stop);
RM_EL = Pm2_EL/Wg2_EL*pi/180
tr0100 = (pi-acos(z))/wa
% saveas(gcf,'marges_A_EL.png')
dx = 0.001
ramp = [0:dx:5];
para =  0.5*ramp'.^2;
y_para = lsim(FTBF_EL_PI,para,ramp); % valeur en reponse a la rampe
y_para_diff = para-y_para; % difference avec la rampe
t2_para_EL_A = dx * find(y_para_diff>y_para_diff(end)*0.98); 
    % trouver le point ou la diff est 98% de l'erreur en regime permanent
    % a la parabole
    
y_obj = y_para_diff(end)- y_para_diff(end)*0.02;

figure()
hold on
line([0 5], [y_para_diff(end)*1.02 y_para_diff(end)*1.02],'LineStyle','--');
line([0 5],[y_para_diff(end)*0.98 y_para_diff(end)*0.98],'LineStyle','--');
plot(ramp',y_para_diff)
saveas(gcf,'ramp_EL_A.png')
    


figure()
hold on
plot(ramp',y_para)
plot(ramp',para)

saveas(gcf,'para_A_EL.png')

figure()
hold on
bode(FTBO_EL)
bode(FTBO_EL*AvPh_EL*PI_EL*band_stop)
legend('originale','finale')
saveas(gcf,'HF_minimal_EL.png')
%% verification de la trajectoire
figure()
lsim(FTBF_EL_PI,utrk,ttrk)
rep_traj = lsim(FTBF_EL_PI,utrk,ttrk);

% calcul de la correlation
R_2 = sum((utrk - mean(rep_traj)).^2) / sum((rep_traj - mean(rep_traj)).^2)
saveas(gcf,'Verif_trajectoire_A_EL.png')





