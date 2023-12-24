%{
*  FILE DESCRIPTION
*  -------------------------------------------------------------------------------------------------------------------
*  File:         Compensator_Design.m
*
*  Description:  Control Engineering- Digital Control Systems Using MATLAB to design a compensator
*
*  -------------------------------------------------------------------------------------------------------------------
*  Author:       Omar Amr
*  Date:         21 Dec 2023
%}

%% Initializing 
clc;
clear;
close all;
syms s;

%% Original System (Without Gain Or Compensator)%%
%%% Plant Continious transfer function (in S-domain)
G_p_Numerator = [2, 1];
G_p_Denominator = [0.2 1.2 1 0];
G_p_S = tf(G_p_Numerator, G_p_Denominator);

%%% Dicretize the plant transfer function (Z-domain)
T_s = 0.1;                      % sampling time = 0.1 seconds
G_p_Z = c2d(G_p_S, T_s, 'zoh'); % discretize using zero order hold

%%% Closed Loop transfer function (in Z-domain)
CL_TF_Z = G_p_Z/(1+G_p_Z);

%%% Palnt transfer function in W domain
G_p_W = d2c(G_p_Z, 'Tustin');

%%% Closed Loop transfer function in W domain
CL_TF_W = d2c(CL_TF_Z, 'Tustin');

%%% Gain, Phase margins and Velocity error of the System 
[G_p_w_Gm, G_p_w_Pm, ~,~] = margin(G_p_W); % Gain and Phase Margins
G_p_w_Gm = 20*log10(G_p_w_Gm);
% Velocity error Calculation
PTF = (2*s + 1)/(0.2*s^3 + 1.2*s^2 + s);
K_v = limit(s*PTF, s, 0);                  

%% Outputs and plots (Before Adding Compensator)
figure(1);
margin(G_p_W);
% title('Gain and Phase Margins No Gain No Compensator)');
figure(2);
step(CL_TF_Z,7);
hold on;
step(CL_TF_W,7);
title('Step responce of the System (No Gain No Compensator)');
fprintf('========== Before Applying Compensator Or Gain ==========\n\n');
fprintf('Gain Margin    = %f dB\n', G_p_w_Gm);
fprintf('Phase Margin   = %f degree\n', G_p_w_Pm);
fprintf('Velocity Error = %d \n', K_v);

%% Adding Gain To The System %%
%%% Open Loop Continious transfer function (in S-domain)
K = 10;  % gain
G_p_S_gain = K * G_p_S;

%%% Dicretize Open Loop transfer function (Z-domain)
G_p_Z_gain = c2d(G_p_S_gain, T_s, 'zoh');

%%% Closed Loop transfer function (in Z-domain)
CL_TF_Z_gain = G_p_Z_gain/(1+G_p_Z_gain);

%%% Open Loop transfer function in W domain
G_p_W_gain = d2c(G_p_Z_gain, 'Tustin');

%%% Closed Loop transfer function in W domain
CL_TF_W_gain = d2c(CL_TF_Z_gain, 'Tustin');

%%% Gain, Phase margins and Velocity error of the System 
[G_p_w_gain_Gm, G_p_w_gain_Pm, ~,~] = margin(G_p_W_gain); % Gain and Phase Margins
G_p_w_gain_Gm = 20*log10(G_p_w_gain_Gm);
% Velocity error Calculation
PTF_gain = K * PTF;
K_v_gain = limit(s*PTF_gain, s, 0);        

%% Outputs and plots (After Adding Gain)
figure(3);
margin(G_p_W_gain);
% title('Gain and Phase Margins After Adding Gain, No Compensator');
figure(4);
step(CL_TF_Z_gain,7);
hold on;
step(CL_TF_W_gain,7);
title('Step responce of the System (With Gain, No Compensator)');

fprintf('\n========== After Applying Gain ==========\n\n');
fprintf('Gain Margin    = %f dB\n', G_p_w_gain_Gm);
fprintf('Phase Margin   = %f degree\n', G_p_w_gain_Pm);
fprintf('Velocity Error = %d \n', K_v_gain);
%% Adding Gain & Compensator To The System %%
% sisotool(G_p_W_gain);   % designing compensator
%%% compensator transfer function in W domain
GD_W_Numerator = [2.5, 1];
GD_W_Denominator = [13, 1];
GD_W = tf(GD_W_Numerator, GD_W_Denominator);

% compensator transfer function in Z domain
GD_Z = c2d(GD_W, T_s, 'Tustin');

%%% Dicretized Open Loop transfer function (Z-domain)
G_p_Z_gain_comp = GD_Z * G_p_Z_gain;

%%% Closed Loop transfer function in Z domain
CL_TL_Z_gain_comp = G_p_Z_gain_comp/(1+G_p_Z_gain_comp);

%%% Open Loop transfer function in W domain
G_p_W_gain_comp = d2c(G_p_Z_gain_comp, 'Tustin');


%%% Closed Loop transfer function in W domain
CL_TL_W_gain_comp = d2c(CL_TL_Z_gain_comp, 'Tustin');


%%% Gain, Phase margins
[G_p_W_gain_comp_Gm, G_p_W_gain_comp_Pm, ~,~] = margin(G_p_W_gain_comp); % Gain and Phase Margins
G_p_W_gain_comp_Gm = 20*log10(G_p_W_gain_comp_Gm);

%% Outputs and plots (After Adding Gain)
figure(5);
margin(G_p_W_gain_comp);
% title('Gain and Phase Margins After Adding Gain & Compensator');
figure(6);
step(CL_TL_Z_gain_comp,7);
hold on;
step(CL_TL_W_gain_comp,7);
title('Step responce of the System (With Gain & Compensator)');

fprintf('\n========== After Applying Gain & Compensator ==========\n\n');
fprintf('Gain Margin    = %f dB\n', G_p_W_gain_comp_Gm);
fprintf('Phase Margin   = %f degree\n', G_p_W_gain_comp_Pm);
% fprintf('Velocity Error = %d \n', K_v_gain);
