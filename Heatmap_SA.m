%% ####################################################################################################################
% Code for the paper:
% Optimized Integration of Solar and Battery Systems in Water Distribution Networks
% By Bhatraj Anudeep, PhD;  Elad Solomons, PhD; Mashor Housh, PhD;
% University of Haifa, abhatraj@campus.haifa.ac.il;elad.salomons@gmail.com; mhoush@univ.haifa.ac.il;
%% ####################################################################################################################
% This code requires:
% Developed under Matlab R2023a
%%
% Plotting heatmap of the SA results.	
clc
clear all
close all

%%
fprintf('\n###################################################################################################\n')
fprintf('######### Optimized Integration of Solar and Battery Systems in Water Distribution Networks ####### \n')
fprintf('###################################################################################################\n\n\n')


load SA_Results.mat


for i = 1:11
    % Convert each inner cell array to a numeric row using cell2mat
    obj_Mat(i, :) = cell2mat(obj_cell{i});
    obj_w_Mat(i, :) = cell2mat(obj_w_cell{i});
    operatingcost_Mat(i, :) = cell2mat(operatingcost_cell{i});
    capitialcost_Mat(i, :) = cell2mat(capitalcost_cell{i});
       
end
 a_value = 0:1000:10000; %A_max_cell
 b_value = 0.5:0.1:1.5;     %Q_demand factor
 highlight_x = 5000;
 highlight_y = 1.0;

figure 
heatmap(a_value,b_value,obj_Mat'/1000,GridVisible="off", ColorbarVisible="off", FontSize=12)
colormap(jet);
colorbar;
% Add labels and title
title('Cost of Total System ($)');
xlabel('Maximum area for solar plant (m^2)');
ylabel("Demand factor");

figure
heatmap(a_value,b_value,value(obj_Mat - obj_w_Mat)'/1000,GridVisible="off", ColorbarVisible="off", FontSize=12)
colormap(jet);
colorbar;
% Add labels and title
title('Cost of Power Distribution System ($)');
xlabel('Maximum area for solar plant (m^2)');
ylabel("Demand factor");

figure
heatmap(a_value,b_value,obj_w_Mat'/1000,GridVisible="off", ColorbarVisible="off", FontSize=12)
colormap(jet);
colorbar;
% Add labels and title
title('Cost of Water Distribution System ($)');
xlabel('Maximum area for solar plant (m^2)');
ylabel("Demand factor");

figure
heatmap(a_value,b_value,operatingcost_Mat'/1000,GridVisible="off", ColorbarVisible="off", FontSize=12)
colormap(jet);
colorbar;
% Add labels and title
title('Operating Cost ($)');
xlabel('Maximum area for solar plant (m^2)');
ylabel("Demand factor");

figure
heatmap(a_value,b_value,capitialcost_Mat'/1000,GridVisible="off", ColorbarVisible="off", FontSize=12)
colormap(jet);
colorbar;
% Add labels and title
title('Capital Cost ($)');
xlabel('Maximum area for solar plant (m^2)');
ylabel("Demand factor");

