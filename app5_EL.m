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
marge = 18%20.5;%19;5
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
para = 0.5 * ramp'.^2;
y_para = lsim(FTBF_EL_PI,para,ramp);
y_para_diff = y_para-para;
y_obj = y_para_diff(end)- y_para_diff(end)*0.02

figure()
hold on
line([0 20],[y_obj y_obj])
plot(ramp',y_para_diff)

%% verification de la trajectoire
% lsim(FTBF_EL_PI,utrk,ttrk)




%% Télescope B
figure
margin(FTBO_EL)

K_AvPh_B = 1/abs(evalfr(FTBO_EL, j*wg_B))


[Gm_AZ Pm_AZ wgm wpm] = margin(K_AvPh_B*FTBO_EL)

phi_AvPh_AZ_B = des_PM_B - Pm_AZ;
alpha_AZ_B = ( 1-sind(phi_AvPh_AZ_B) )/( 1+sind(phi_AvPh_AZ_B) );
T_AZ_B = 1/(wg_B*sqrt(alpha_AZ_B));

z = -1/T_AZ_B;
p = -1/(T_AZ_B*alpha_AZ_B);

Ka_AvPh_B = K_AvPh_B/sqrt(alpha_AZ_B);
s = tf('s');
G_AvPh_B = Ka_AvPh_B*(s-z)/(s-p);

FTBO_EL_B1 = series(G_AvPh_B, FTBO_EL);
[num den] = tfdata(FTBO_EL_B1, 'v');
Kvel = num(end)/den(end-1);
erp_ramp_AZ_B = 1/Kvel

% figure
% margin(FTBO_EL_B1)

K_RePh_B = (1/des_erp_B)/Kvel;

% figure
% margin(K_RePh_B*FTBO_EL_B1)

T_AZ_Re_B = 10/wg_B

z = -1/T_AZ_Re_B
z = -0.68
p = -1/(K_RePh_B*T_AZ_Re_B)

Kr_AZ_B = 1/abs((j*wg_B - z)/(j*wg_B-p))
Kr_AZ_B = 0.97

G_RePh_B = Kr_AZ_B*(s-z)/(s-p)
FTBO_EL_B2 = series(G_RePh_B, FTBO_EL_B1)

[num den] = tfdata(FTBO_EL_B2, 'v');
Kvel = num(end)/den(end-1);
erp_ramp_AZ_B = 1/Kvel

figure
margin(FTBO_EL_B2)



figure
step(feedback(FTBO_EL_B2,1), [0:0.001:5])

w_c = 123 %pic
w_width = 10
num_band_stop = [1 0 w_c^2];
den_band_stop = [1 w_width w_c^2];
band_stop = tf(num_band_stop, den_band_stop)

FTBO_EL_B3 = series(FTBO_EL_B2, band_stop)

[Gm_AZ Pm_AZ wgm wpm] = margin(FTBO_EL_B3);
zeta = sqrt(tand(Pm_AZ)*sind(Pm_AZ))/2;
wn = wpm*tand(Pm_AZ)/(2*zeta);
ts = 4/(zeta*wn)

BW = wpm*sqrt(1-zeta^2+sqrt(4*zeta^4-4*zeta^2+2))/sqrt(sqrt(1+4*zeta^4)-2*zeta^2)

figure
margin(FTBO_EL_B3)
%xlim([40 60])

figure
step(feedback(FTBO_EL_B3,1), [0:0.001:5])

figure
rlocus(FTBO_EL_B3)

t = [0:0.01:30];
u = t;
figure
lsim(feedback(FTBO_EL_B3, 1),u,t) 