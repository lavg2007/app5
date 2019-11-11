%traduction des specs
theta = atan(-pi/log(des_MP_A/100))
z = cos(theta)
wn_1 = (4/des_ts2_A)/z;
wn_2 = (1+1.1*z+1.4*z^2)/des_tm1090_A;
wn_3 = (pi-acos(z))/(des_tm0100_A * sqrt(1-z^2));
wn = max([wn_1 wn_2 wn_3]);
wa = wn*sin(theta);
% poles desires
s(1) = -z*wn + wa*i;
s(2) = -z*wn - wa*i;

%% B

z_B = sqrt(tand(des_PM_B)*sind(des_PM_B))/2
wg_B = des_BW_B*sqrt(sqrt(1+4*z_B^4)-2*z_B^2)/sqrt((1-2*z_B^2)+sqrt(4*z_B^4-4*z_B^2+2))