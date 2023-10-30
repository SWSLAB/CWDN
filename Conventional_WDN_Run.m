%% ####################################################################################################################
% Code for the paper:
% Optimized Integration of Solar and Battery Systems in Water Distribution Networks
% By Bhatraj Anudeep, PhD; Mashor Housh, PhD; Elad Salomons, PhD
% University of Haifa, abhatraj@campus.haifa.ac.il;mhoush@univ.haifa.ac.il;elad.salomons@gmail.com
%% ####################################################################################################################
% This code requires:
% YALMIP toolbox: https://yalmip.github.io/
% Optitool box: https://inverseproblem.co.nz/OPTI/
% Gurobi:https://www.gurobi.com/features/academic-named-user-license/
% Developed under Matlab R2023a

% Determine the least cost design the conventional WDN system.	
clc
clear all
close all
yalmip('clear')

%%
fprintf('\n###################################################################################################\n')
fprintf('######### Optimized Integration of Solar and Battery Systems in Water Distribution Networks ####### \n')
fprintf('###################################################################################################\n\n\n')

%% Define parameters

Define_Parameters
eta_c = 0;	                % Efficiency of battery on charging .
eta_d = 0;	                % Efficiency of battery on discharging.
Amax = 0;                   % Maximum area available for the solar panels (m^2).

%% Formulate and solve the model
CWDN_Model

%% Summarize and plot results
Summarize_and_Plot_Results


 