clear all;
close all;
clc;
%% trouver les caracteristiques du systeme
kp = 1
num = 1
den = [1 3 2]
g = tf(num,den)
FTBF.FT = feedback(g,kp)
[FTBF.num, FTBF.den] = tfdata(FTBF.FT,'v')
FTBF.wn = sqrt(FTBF.den(3))
FTBF.z = FTBF.den(2)/(2*FTBF.wn)

P_ori(1) = -FTBF.z*FTBF.wn + FTBF.wn*sqrt(1-FTBF.z^2)*i;
P_ori(2) = -FTBF.z*FTBF.wn - FTBF.wn*sqrt(1-FTBF.z^2)*i;
Kpos_voulue = kp /((0+1)*(0+2));
Kpos = 1/0.2 - 1;
Kpos_prime = Kpos/Kpos_voulue;
i = 0
T = [0:0.1:30];
U = ones(size(T));
% y(end)>0.9 && y(end)<1.1)
for d = 1:0.01:5;
i = i+1;
z = real(P_ori(1))/(d);
p = z/Kpos_prime;

% creation du compensateur
s = tf('s');
gc = (s-z)/(s-p);

y(i,:) = lsim(feedback(gc*g,kp),U,T);
end
% [v,idx] = min(abs(err))

% figure(2)
% plot(y(idx,:))
%% num 1 lab
clear all;
close all;
clc;
tr = 0.004
tp = 0.008
num = [4500]
den = [1 361.2 0]
FBTO = tf(num,den)
% traduire en poles desires
theta = atan(-pi/log(6/100))
z = cos(theta)
wn_1 = (4/0.010)/z'
wn_2 = (1+1.1*z+1.4*z^2)/tr;
wn_3 = pi/(tp*sqrt(1-z^2));
wn = wn_1;
P_des(1) = -z*wn + wn*sqrt(1-z^2)*i;
P_des(2) = -z*wn - wn*sqrt(1-z^2)*i;

% creation d'un avance de phase
phase_G = rad2deg(angle(num/(P_des(1)*(P_des(1)+361.2))));
delta_phi = -180 - phase_G+360;
phi = 180 - rad2deg(atan2(imag(P_des(1)),real(P_des(1))));
alpha = 180-phi;
phi_z = (alpha + delta_phi)/2;
phi_p = (alpha - delta_phi)/2;
z = real(P_des(1)) - imag(P_des(1))/tan(deg2rad(phi_z));
p = real(P_des(1)) - imag(P_des(1))/tan(deg2rad(phi_p));

s = P_des(1)
ka = 1/norm((s-z)/(s-p)*num/(s*(s+361)));
s = tf('s');
AvPh = ka*(s-z)/(s-p);

% creation du retard de phase
H = AvPh*FBTO;
H = minreal(H*s);
[H_num,H_den] = tfdata(H,'v');
k_vel_act = (H_num(end)/H_den(end));
k_vel_des = 1/0.00005;
k_vel = k_vel_des/k_vel_act;
z_pi = real(P_des(1))/10;
p_pi = z_pi/k_vel;
s = P_des(1)
kr = 1 % comme ecrit dans les notes
s = tf('s');
RePh = kr*(s-z_pi)/(s-p_pi);

step(feedback(RePh*AvPh*FBTO,1))
FTBF = minreal(feedback(RePh*AvPh*FBTO,1))

% creation du PD
