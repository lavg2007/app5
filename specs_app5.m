%% Script de l'app 5 S5e (section du materiel fourni)
% Par Hubert Dube
% Debute le 7/11/2019
%% Fonctions de transfert des deux telescope
% Azimut
numAZ = [1.59e09]
denAZ = [1 1020.51 25082.705 3102480.725 64155612.5 82700000 0]
FTBO_AZ = tf(numAZ,denAZ)

% Elevation
numEL = [7.95e09]
denEL = [1 1020.51 37082.705 15346520.725 320776412.5 413500000 0]
FTBO_EL = tf(numEL,denEL)
%% Specification du telescope A
% critere de design
des_MP_A = 25%;13.8 % faire ca dans un FOR 
des_ts2_A = 1;
des_tm1090_A = 0.15;
des_tm0100_A = 0.25;
des_erp_A = 0;       % Erreur en Regime Permanent
des_eru_AZ_A = 0.05; % Erreur Rampe Unitaire en Regime Permanent
des_eru_EL_A = 0;
des_epu_AZ_A = NaN;  % Erreur Parabole Unitaire en Regime Permanent
des_epu_EL_A = 0.1;

% critere de securite
sec_PM_A = 10;              % db
sec_RM_AZ_A = 0.10;         % s
sec_RM_EL_A = 0.08;         % s  
sec_Atten_Vib_AZ_A = -15;   % db
sec_Atten_Vib_AZ_A = -15;   % db

% critere d'acceptation finale
acc_MP_A = 30;
acc_ts2_eu_AZ_A = 1.25; % Erreur en Regime Permanent Echellon unitaire
acc_ts2_eu_EL_A = 1.50;
acc_ts2_ru_AZ_A = 1.25; % Erreur Rampe Unitaire en Regime Permanent
acc_ts2_ru_EL_A = 1.50;
acc_ts2_pu_AZ_A = NaN;  % Erreur Parabole Unitaire en Regime Permanent
acc_ts2_pu_EL_A = 3;
%% Specification du telescope B
% critere de design
des_BW_B = 10;
des_PM_B = 50; % +- 5 deg
des_erp_B = 0.005;
sec_Atten_Vib_B = -12;   % db
% sec_14db_B = TDB 

% critere de securite
sec_GM = 15; %db

% critere d'acceptation finale
acc_ts2_ru_B = 14; %s
