%% Telescope B
specs_app5
trad_specs
Profile_Tracking;


figure
margin(FTBO_AZ)

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
% z = -0.68
p = -1/(K_RePh_B*T_AZ_Re_B)

Kr_AZ_B = 1/abs((j*wg_B - z)/(j*wg_B-p));
Kr_AZ_B = 1;


G_RePh_B = Kr_AZ_B*(s-z)/(s-p)
FTBO_AZ_B2 = series(G_RePh_B, FTBO_AZ_B1)



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

[num den] = tfdata(FTBO_AZ_B3, 'v');
Kvel = num(end)/den(end-1);
erp_ramp_AZ_B = 1/Kvel

[Gm_AZ Pm_AZ wgm wpm] = margin(FTBO_AZ_B3);
zeta = sqrt(tand(Pm_AZ)*sind(Pm_AZ))/2;
wn = wpm*tand(Pm_AZ)/(2*zeta);
ts = 4/(zeta*wn)
BW = wpm*sqrt(1-zeta^2+sqrt(4*zeta^4-4*zeta^2+2))/sqrt(sqrt(1+4*zeta^4)-2*zeta^2)

figure
margin(FTBO_AZ_B3)
saveas(gcf,[pwd '\Rapport\figures\B_AZ_margin_bandstop.png'])

figure
step(feedback(FTBO_AZ_B3,1), [0:0.001:5])

figure
rlocus(FTBO_AZ_B3)

%%

t = [0:0.01:30];
u = t;

y3 = lsim(feedback(FTBO_AZ_B3, 1),u,t); 

figure
plot(t,y3'-u)
saveas(gcf,[pwd '\Rapport\figures\end_erreur_B_AZ.png'])

FTBF_AZ_B = feedback(FTBO_AZ_B3,1)
%% verification de la trajectoire

figure()
lsim(FTBF_AZ_B,utrk,ttrk)
rep_traj_B = lsim(FTBF_AZ_B,utrk,ttrk);

% calcul de la correlation
R_2_B = sum((utrk - mean(rep_traj_B)).^2) / sum((rep_traj_B - mean(rep_traj_B)).^2)
saveas(gcf,[pwd '\Rapport\figures\Verif_trajectoire_B_AZ.png'])