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

mrg_PM_AZ = 5
% mrg_PM_AZ = 5

mrg_PM_EL = 5
% mrg_PM_EL = 5

z_B_AZ = sqrt(tand(des_PM_B+mrg_PM_AZ)*sind(des_PM_B+mrg_PM_AZ))/2
z_B_EL = sqrt(tand(des_PM_B+mrg_PM_EL)*sind(des_PM_B+mrg_PM_EL))/2
wg_B_AZ = des_BW_B*sqrt(sqrt(1+4*z_B_AZ^4)-2*z_B_AZ^2)/sqrt((1-2*z_B_AZ^2)+sqrt(4*z_B_AZ^4-4*z_B_AZ^2+2))
wg_B_EL = des_BW_B*sqrt(sqrt(1+4*z_B_EL^4)-2*z_B_EL^2)/sqrt((1-2*z_B_EL^2)+sqrt(4*z_B_EL^4-4*z_B_EL^2+2))