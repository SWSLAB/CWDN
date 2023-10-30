
%% ####################################################################################################################
% Code for the paper:
% Optimized Integration of Solar and Battery Systems in Water Distribution Networks
% By Bhatraj Anudeep, PhD;  Elad Salomons, PhD; Mashor Housh, PhD;
% University of Haifa, abhatraj@campus.haifa.ac.il;elad.salomons@gmail.com;mhoush@univ.haifa.ac.il
%% ####################################################################################################################
% This code requires:
% YALMIP toolbox: https://yalmip.github.io/
% Optitool box: https://inverseproblem.co.nz/OPTI/
% Gurobi:https://www.gurobi.com/features/academic-named-user-license/
% Developed under Matlab R2023a

% Run the CWDN model for Sensitivity Analysis
clc
clear all
close all
yalmip('clear')
%%
fprintf('\n###################################################################################################\n')
fprintf('######### Optimized Integration of Solar and Battery Systems in Water Distribution Networks ####### \n')
fprintf('###################################################################################################\n\n\n')

%% Define paraemters
Define_Parameters
Q_dem_act = Q_dem;
Amax = 0:1000:10000;         % Loop for maximum area available for the solar panels (m^2).
factor = 0.5:0.1:1.5;           % Loop for demand factor.
obj_w_cell = cell(1,0);
obj_p_cell = cell(1,0);
obj_cell = cell(1,0);
capitalcost_cell = cell(1,0);
operatingcost_cell = cell(1,0);

counter = 0;
numIt_Outer = numel(Amax);
numIt_Inner = numel(factor);
totalIterations = numIt_Outer * numIt_Inner;
i = 1:11;
for Amax = 0:1000:10000
    % Initialize a cell array to store the values of E_BST_max for each factor
    obj_w_cell_Amax = cell(1,0);
    obj_p_cell_Amax = cell(1,0);
    obj_cell_Amax = cell(1,0);
    capitalcost_Amax = cell(1,0);
    operatingcost_Amax = cell(1,0);
    for  factor = 0.5:0.1:1.5
        Q_dem= factor*Q_dem_act;
        
        % Formulate and solve the model
        CWDN_Model
        
        % Save values in cell arrays
        obj_w_cell_Amax{end+1} = value (obj_W);
        obj_p_cell_Amax{end+1} = value (obj_P);
        obj_cell_Amax{end+1} = value(obj);
        capitalcost_Amax{end+1} = value (obj_W + COS(A_p*N_pan,P_S_std) + COB(E_BST_max));
        operatingcost_Amax{end+1} = value (obj - (obj_W + COS(A_p*N_pan,P_S_std) + COB(E_BST_max)));
        
        counter = counter + 1;
        disp(['Iteration: ' num2str(counter) ' out of ' num2str(totalIterations)]);
    end
    obj_w_cell{end+1} = obj_w_cell_Amax;
    obj_p_cell{end+1} = obj_p_cell_Amax;
    obj_cell{end+1} = obj_cell_Amax;
    capitalcost_cell{end+1} = capitalcost_Amax;
    operatingcost_cell{end+1} = operatingcost_Amax;
end

% Save the E_BST_max_cell to a mat file
save('SA_Results.mat', "capitalcost_cell", "operatingcost_cell", "obj_cell", 'obj_p_cell', 'obj_w_cell');

