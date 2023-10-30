%% ####################################################################################################################
% Code for the paper:
% Optimized Integration of Solar and Battery Systems in Water Distribution Networks
% By Bhatraj Anudeep, PhD; Mashor Housh, PhD; Elad Salomons, PhD
% University of Haifa, abhatraj@campus.haifa.ac.il;mhoush@univ.haifa.ac.il;elad.salomons@gmail.com
%% ####################################################################################################################
% This code requires:
% Developed under Matlab R2023a

%% General  Parameters:

T= 24;                  % Hours in a day.
SEA = 4;                % No. of seasons in a year.
W = 2;                  % Weekend or weekdays.

%% data files

load ('rate_tariff_s_w.mat')        % Energy per unit cost (Buying Price)(NIS/kWh).
E_tariff = (rate_tariff/3.5)*1.17 ; % NIS to $ and add VAT of 17%.
load ('n_s_w.mat');                 % Number of days in every season and week types.
tmp(1,:,:) = n_s_w;
n_s_w_mat = repmat(tmp,T,1,1);
load ('demand_s_w.mat');
Q_dem = demand_s_w;             % Demand time series (m^3/hr).
load('irradiation_s_w.mat');
IR_t = irradiation/1000;            % Converting into (kW/m^2);

%% Water parameters

gama = 9810;       % Specific weight of water (N/m^3).
C = 130;           % Hazen-Williams Coefficient.
l= 1450;           % Length of the pipe (m).
H_tank = 275;	   % Head of tank (m).
H_res =55 ;        % Head of reservoir (m).
d_min = 0.1;       % Minimum diameter of the pipe (m).
d_max= 0.7	;      % Maximum diameter of the pipe (m).
eta = 0.85;        % Efficiency of the pump.
dt = 1;            % Time step duration (hr).
r_w = 0.055;       % Rate of interest for WDN (-).
N_pipe = 50;       % Life span of the pipes (yrs).
N_tank = 50;       % Life span of the tank (yrs).
N_pump = 15;       % Life span of the pump (yrs).

V_min= [300, 300; 150, 150; 150, 150; 150, 150] ;      % Minimum volume of the tanks (m^3).

% Capital recovery factor https://en.wikipedia.org/wiki/Capital_recovery_factor;
Crf_pipe =  (r_w *(1 +r_w)^N_pipe)/ ((1 +r_w)^N_pipe -1);
Crf_tank =  (r_w *(1 +r_w)^N_tank)/ ((1 +r_w)^N_tank -1);
Crf_pump =  (r_w *(1 +r_w)^N_pump)/ ((1 +r_w)^N_pump -1);

%% Grid parameters

E_S_tariff=  E_tariff/2;        % Selling price of energy defined as half of the buying price($/kWh);
eta_c = 0.95;	                % Efficiency of battery on charging.
eta_d = 0.95;	                % Efficiency of battery on discharging.
SOC_max = 0.9;                  % Maximum State of charge.
SOC_min = 0.1;                  % Minimum State of charge.
DF = 0.8;                       % Derating factor of solar panels.
IR_std = 1;                     % Standard irradiance (kW/m^2).
eta_pan = 0.18;	                % Efficiency of panel.
eta_inv = 0.95;                 % Efficiency of the inverter.
A_p = 1.8;                      % Area of the solar panel Standard (m^2).
N_bat = 20;                     % Life span of the battery storage (yrs).
N_sol = 20;                     % Life span of solar panels (yrs).
r_p = r_w;                      % Rate of interest for PDN.
C_c_time = 2;                   % Maximum charging time of battery source(hr).
C_d_time = 4;                   % Maximum discharging time of battery source(hr).
Amax = 5000;                    % Maximum area available for the solar panels (m^2).

% Capital recovery factor https://en.wikipedia.org/wiki/Capital_recovery_factor;
Crf_solar = (r_p/(1-(1+r_p)^-N_sol));
Crf_bat =   (r_p/(1-(1+r_p)^-N_bat));

