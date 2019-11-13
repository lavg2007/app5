close all
clear all
clc
%% Script de l'app 5 S5e
% Par Hubert Dube
% Debute le 7/11/2019
specs_app5
%% Telescope A ------------------------------------------------------------
%--------------------------------------------------------------------------
%-------------------------------AZIMUT-------------------------------------
%--------------------------------------------------------------------------
%% traduction des specifications
trad_specs
Profile_Tracking;

% creation du lieu des raciines
% figure(1)
% hold on
% rlocus(FTBO_AZ)
% scatter(real(s),imag(s),'p')
% title('LdR FTBO-AZ')
% reponse a l'echellon initiale
figure()
step(feedback(FTBO_AZ,1))
saveas(gcf,'FTBO_AZ.png')
figure()
margin(FTBO_AZ)
saveas(gcf,'marges_FTBO_AZ.png')
%% Compensateur pour AZIMUT
% creation d'un avance de phase
phase_AZ = rad2deg(angle(numAZ/polyval(denAZ,s(1))))
delta_phi_AvPh_AZ = -180 - phase_AZ + 360
phi_AvPh_AZ = 180 - rad2deg(atan2(imag(s(1)),real(s(1))));
alpha_AvPh_AZ = 180-phi_AvPh_AZ;
phi_z_AvPh_AZ = (alpha_AvPh_AZ + delta_phi_AvPh_AZ)/2;
phi_p_AvPh_AZ = (alpha_AvPh_AZ - delta_phi_AvPh_AZ)/2;
z_AvPh_AZ = real(s(1)) - imag(s(1))/tan(deg2rad(phi_z_AvPh_AZ));
p_AvPh_AZ = real(s(1)) - imag(s(1))/tan(deg2rad(phi_p_AvPh_AZ));


ka_AvPh_AZ = 1/norm((s(1)-z_AvPh_AZ)/(s(1)-p_AvPh_AZ)* numAZ/polyval(denAZ,s(1)))
num_AvPh_AZ = [1 -z_AvPh_AZ];
den_AvPh_AZ = [1 -p_AvPh_AZ];
AvPh_AZ = ka_AvPh_AZ*tf(num_AvPh_AZ,den_AvPh_AZ);
% ka_AvPh2_AZ = 1/norm(ka_AvPh_AZ*(s(1)-z_AvPh_AZ)/(s(1)-p_AvPh_AZ)*(s(1)-z_AvPh_AZ)/(s(1)-p_AvPh_AZ)* numAZ/polyval(denAZ,s(1)))
% AvPh2_AZ = ka_AvPh2_AZ*tf(num_AvPh_AZ,den_AvPh_AZ);
% creation du lieu des racines avec l'AvPH
figure(2)
hold on
rlocus(FTBO_AZ*AvPh_AZ)
scatter(real(s),imag(s),'p')
title('LdR FTBO-AZ*AvPh-AZ')


figure(3)
hold on
step(feedback(FTBO_AZ*AvPh_AZ,1),5)
step(feedback(FTBO_AZ,1),5)
[y_FTBF_AvPh_AZ,t] = step(feedback(FTBO_AZ*AvPh_AZ,1),5);
y_FTBF_AZ = step(feedback(FTBO_AZ,1),5);
legend('FTBO-AZ*AvPh','FTBO-AZ');
title('Step comparaison FTBO-AZ*AvPh-AZ');

figure(30)
hold on
margin(FTBO_AZ*AvPh_AZ)
saveas(gcf,'FTBO_AZ_AvPh_marge.png')
%% verification du systeme asservi par l'AvPh
FTBF_AZ_AvPh = feedback(FTBO_AZ*AvPh_AZ,1);
[num_FTBO_AZ_AvPh,den_FTBO_AZ_AvPh] = tfdata(FTBO_AZ*AvPh_AZ,'v');
stepinfo(FTBF_AZ_AvPh)
marge = 2
% reponse a une rampe unitaire
Kvel = polyval(num_FTBO_AZ_AvPh,0)/polyval([den_FTBO_AZ_AvPh(1:end-1)],0);
eru_AZ_A = 1/Kvel

% verification des marges de securite
figure(5)
margin(FTBO_AZ*AvPh_AZ)
[Gm_AZ,Pm_AZ,Wp_AZ,Wg_AZ] = margin(FTBO_AZ*AvPh_AZ);
RM_AZ = Pm_AZ/Wg_AZ*pi/180


% creation d'un second AvPh pour respecter la marge de retard
PM_AZ_des = sec_RM_AZ_A *180/pi *Wg_AZ
K_AvPh2 = 1/norm(polyval(num_FTBO_AZ_AvPh,0+Wg_AZ*i)/polyval([den_FTBO_AZ_AvPh],0+Wg_AZ*i))

pm = angle(K_AvPh2*polyval(num_FTBO_AZ_AvPh,0+Wg_AZ*i)/polyval([den_FTBO_AZ_AvPh],0+Wg_AZ*i))*180/pi - -180;
delta_phi_AvPh2 = PM_AZ_des - pm + marge

alpha_AvPh2 = (1-sind(delta_phi_AvPh2)) / (1+sind(delta_phi_AvPh2))
T_AvPh2 = 1/(Wg_AZ*sqrt(alpha_AvPh2))
z_AvPh2 = -1/T_AvPh2
p_AvPh2 = -1/(alpha_AvPh2*T_AvPh2)
ka_AvPh2 = K_AvPh2/ sqrt(alpha_AvPh2)

num_AvPh2_AZ = [1 -z_AvPh2];
den_AvPh2_AZ = [1 -p_AvPh2];
AvPh2_AZ = ka_AvPh2*tf(num_AvPh2_AZ,den_AvPh2_AZ);

%% verification avec nouvelle compensation
figure(6)
margin(FTBO_AZ*AvPh_AZ*AvPh2_AZ)
[Gm2_AZ,Pm2_AZ,Wp2_AZ,Wg2_AZ] = margin(FTBO_AZ*AvPh_AZ*AvPh2_AZ);
RM2_AZ = Pm2_AZ/Wg2_AZ*pi/180

figure(7)
hold on
step(feedback(FTBO_AZ,1),5)
step(feedback(FTBO_AZ*AvPh_AZ,1),5)
step(feedback(FTBO_AZ*AvPh_AZ*AvPh2_AZ,1),5)
legend('FTBO_AZ','AvPh_AZ','AvPh2_AZ')
[num_FTBO_AZ_AvPh2,den_FTBO_AZ_AvPh2] = tfdata(FTBO_AZ*AvPh_AZ*AvPh2_AZ,'v')
Kvel2 = polyval(num_FTBO_AZ_AvPh2,0)/polyval([den_FTBO_AZ_AvPh2(1:end-1)],0);
eru2_AZ_A = 1/Kvel2

%% ajout d'un coupe bande
% trouver la frequence de coupure
freq_coup = 54.8 % rad/sec trouve avec bode
w_width = 12
num_band_stop = [1 0 freq_coup^2];
den_band_stop = [1 w_width freq_coup^2];
band_stop = tf(num_band_stop,den_band_stop);

figure(9);
hold on;
step(feedback(FTBO_AZ*AvPh_AZ*AvPh2_AZ*band_stop,1),5);
step(feedback(FTBO_AZ*AvPh_AZ*AvPh2_AZ,1),5);
step(feedback(FTBO_AZ*AvPh_AZ,1),5);
step(feedback(FTBO_AZ,1),5);
legend('band_stop','AvPh2','AvPh1','ori');
saveas(gcf,'step_A_AZ_comp_wcb.png')
xlim([0 3])
stepinfo(feedback(FTBO_AZ*AvPh_AZ*AvPh2_AZ*band_stop,1))
%% verification final
figure(10)
margin(FTBO_AZ*AvPh_AZ*AvPh2_AZ*band_stop)
[Gm3_AZ,Pm3_AZ,Wp3_AZ,Wg3_AZ] = margin(FTBO_AZ*AvPh_AZ*AvPh2_AZ*band_stop);
RM3_AZ = Pm3_AZ/Wg3_AZ*pi/180
saveas(gcf,'marges_A_AZ.png')
[num_FTBO_AZ_AvPh3,den_FTBO_AZ_AvPh3] = tfdata(FTBO_AZ*AvPh_AZ*AvPh2_AZ*band_stop,'v');
Kvel3 = polyval(num_FTBO_AZ_AvPh3,0)/polyval([den_FTBO_AZ_AvPh3(1:end-1)],0);
eru3_AZ_A = 1/Kvel3

stepinfo(feedback(FTBO_AZ*AvPh_AZ*AvPh2_AZ*band_stop,1))
FTBF_AZ = feedback(FTBO_AZ*AvPh_AZ*AvPh2_AZ*band_stop,1)
ramp = [0:0.1:20];
y_ramp = lsim(FTBF_AZ,ramp,ramp);
y_ramp_diff = y_ramp-ramp';

figure()
hold on
plot(ramp',ramp')
plot(ramp',y_ramp)
saveas(gcf,'ramp_AZ_A.png')

figure()
hold on
plot(ramp',ramp')
plot(ramp',y_ramp)
xlim([19 20]);
ylim([19 20]);
saveas(gcf,'ramp_AZ_A_zoom.png')

%% verification de la trajectoire
figure()
lsim(FTBF_AZ,utrk,ttrk)
rep_traj_A = lsim(FTBF_AZ,utrk,ttrk);

% calcul de la correlation
R_2_A = sum((utrk - mean(rep_traj_A)).^2) / sum((rep_traj_A - mean(rep_traj_A)).^2)
saveas(gcf,'Verif_trajectoire_A_AZ.png')


%% Telescope B

% figure
% margin(FTBO_AZ)

K_AvPh_B = 1/abs(evalfr(FTBO_AZ, j*wg_B))


[Gm_AZ Pm_AZ wgm wpm] = margin(K_AvPh_B*FTBO_AZ)

phi_AvPh_AZ_B = des_PM_B - Pm_AZ;
alpha_AZ_B = ( 1-sind(phi_AvPh_AZ_B) )/( 1+sind(phi_AvPh_AZ_B) );
T_AZ_B = 1/(wg_B*sqrt(alpha_AZ_B));

z = -1/T_AZ_B;
p = -1/(T_AZ_B*alpha_AZ_B);

Ka_AvPh_B = K_AvPh_B/sqrt(alpha_AZ_B);
s = tf('s');
G_AvPh_B = Ka_AvPh_B*(s-z)/(s-p);

FTBO_AZ_B1 = series(G_AvPh_B, FTBO_AZ);
[num den] = tfdata(FTBO_AZ_B1, 'v');
Kvel = num(end)/den(end-1);
erp_ramp_AZ_B = 1/Kvel

% figure
% margin(FTBO_AZ_B1)

K_RePh_B = (1/des_erp_B)/Kvel;

% figure
% margin(K_RePh_B*FTBO_AZ_B1)

T_AZ_Re_B = 10/wg_B

z = -1/T_AZ_Re_B
z = -0.68
p = -1/(K_RePh_B*T_AZ_Re_B)

Kr_AZ_B = 1/abs((j*wg_B - z)/(j*wg_B-p));
% Kr_AZ_B = 0.95;


G_RePh_B = Kr_AZ_B*(s-z)/(s-p)
FTBO_AZ_B2 = series(G_RePh_B, FTBO_AZ_B1)

[num den] = tfdata(FTBO_AZ_B2, 'v');
Kvel = num(end)/den(end-1);
erp_ramp_AZ_B = 1/Kvel

figure
margin(FTBO_AZ_B2)




figure
step(feedback(FTBO_AZ_B2,1), [0:0.001:5])

w_c = 54.8 %pic
w_width = 10
num_band_stop = [1 0 w_c^2];
den_band_stop = [1 w_width w_c^2];
band_stop = tf(num_band_stop, den_band_stop)
[num_cheb den_cheb] = cheby1(2,0,[w_c-w_width/2 w_c+w_width/2], 'stop', 's')
band_stop = tf(num_cheb, den_cheb)


FTBO_AZ_B3 = series(FTBO_AZ_B2, band_stop)

[Gm_AZ Pm_AZ wgm wpm] = margin(FTBO_AZ_B3);
zeta = sqrt(tand(Pm_AZ)*sind(Pm_AZ))/2;
wn = wpm*tand(Pm_AZ)/(2*zeta);
ts = 4/(zeta*wn)
BW = wpm*sqrt(1-zeta^2+sqrt(4*zeta^4-4*zeta^2+2))/sqrt(sqrt(1+4*zeta^4)-2*zeta^2)

figure
margin(FTBO_AZ_B3)
xlim([40 60])

figure
step(feedback(FTBO_AZ_B3,1), [0:0.001:5])

figure
rlocus(FTBO_AZ_B3)

t = [0:0.01:30];
u = t;
figure
lsim(feedback(FTBO_AZ_B3, 1),u,t) 
FTBF_AZ_B = feedback(FTBO_AZ_B3,1)
%% verification de la trajectoire

figure()
lsim(FTBF_AZ_B,utrk,ttrk)
rep_traj_B = lsim(FTBF_AZ_B,utrk,ttrk);

% calcul de la correlation
R_2_B = sum((utrk - mean(rep_traj_B)).^2) / sum((rep_traj_B - mean(rep_traj_B)).^2)
saveas(gcf,'Verif_trajectoire_B_AZ.png')