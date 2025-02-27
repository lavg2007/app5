%% T�lescope B
specs_app5;
trad_specs;
Profile_Tracking;

%% Analyse du syst�me
figure
margin(FTBO_EL)


%% Conception de l'avance de phase
K_AvPh_B = 1/abs(evalfr(FTBO_EL, j*wg_B))


[Gm_EL Pm_EL wgm wpm] = margin(K_AvPh_B*FTBO_EL)

phi_AvPh_EL_B = des_PM_B - Pm_EL;
alpha_EL_B = ( 1-sind(phi_AvPh_EL_B) )/( 1+sind(phi_AvPh_EL_B) );
T_EL_B = 1/(wg_B*sqrt(alpha_EL_B));

z = -1/T_EL_B;
p = -1/(T_EL_B*alpha_EL_B);

Ka_AvPh_B = K_AvPh_B/sqrt(alpha_EL_B);
s = tf('s');
G_AvPh_B = Ka_AvPh_B*(s-z)/(s-p);

FTBO_EL_B1 = series(G_AvPh_B, FTBO_EL);
[num den] = tfdata(FTBO_EL_B1, 'v');
Kvel = num(end)/den(end-1);
erp_ramp_EL_B = 1/Kvel
% figure
% margin(FTBO_EL_B1)


%% Conception du Retard de phase

K_RePh_B = (1/des_erp_B)/Kvel;

% figure
% margin(K_RePh_B*FTBO_EL_B1)

T_EL_Re_B = 10/wg_B

z = -1/T_EL_Re_B
p = -1/(K_RePh_B*T_EL_Re_B)

% Kr_EL_B = 1/abs((j*wg_B - z)/(j*wg_B-p))
Kr_EL_B = 1

G_RePh_B = Kr_EL_B*(s-z)/(s-p)
FTBO_EL_B2 = series(G_RePh_B, FTBO_EL_B1)

[num den] = tfdata(FTBO_EL_B2, 'v');
Kvel = num(end)/den(end-1);
erp_ramp_EL_B = 1/Kvel

figure
margin(FTBO_EL_B2)

figure
step(feedback(FTBO_EL_B2,1), [0:0.001:5])

%% Coupe Bande
w_c = 123 %pic
w_width = 10
num_band_stop = [1 0 w_c^2];
den_band_stop = [1 w_width w_c^2];
band_stop = tf(num_band_stop, den_band_stop)
[num_cheb den_cheb] = cheby1(2,0,[w_c-w_width/2 w_c+w_width/2], 'stop', 's')
band_stop = tf(num_cheb, den_cheb)


FTBO_EL_B3 = series(FTBO_EL_B2, band_stop)

%% Analyse du design initial

[Gm_EL Pm_EL wgm wpm] = margin(FTBO_EL_B3);
zeta = sqrt(tand(Pm_EL)*sind(Pm_EL))/2;
wn = wpm*tand(Pm_EL)/(2*zeta);
ts = 4/(zeta*wn)

BW = wpm*sqrt(1-zeta^2+sqrt(4*zeta^4-4*zeta^2+2))/sqrt(sqrt(1+4*zeta^4)-2*zeta^2)
figure
margin(FTBO_EL_B3)
saveas(gcf,[pwd '\Rapport\figures\B_EL_margin_bandstop.png'])

%xlim([40 60])

figure
step(feedback(FTBO_EL_B3,1), [0:0.001:5])

figure
rlocus(FTBO_EL_B3)

t = [0:0.01:30];
u = t;
y3 = lsim(feedback(FTBO_EL_B3, 1),u,t) 

figure
plot(t, y3'-u)
saveas(gcf,[pwd '\Rapport\figures\B_EL_error_final.png'])



FTBF_EL_B=  feedback(FTBO_EL_B3,1)


%% verification de la trajectoire
figure()
lsim(FTBF_EL_B,utrk,ttrk)

rep_traj_B = lsim(FTBF_EL_B,utrk,ttrk);

% calcul de la correlation
R_2_B = sum((utrk - mean(rep_traj_B)).^2) / sum((rep_traj_B - mean(rep_traj_B)).^2)
saveas(gcf,'Verif_trajectoire_B_EL.png')