%% ####################################################################################################################
% Code for the paper:
% Optimized Integration of Solar and Battery Systems in Water Distribution Networks
% By Bhatraj Anudeep, PhD; Mashor Housh, PhD; Elad Salomons, PhD
% University of Haifa, abhatraj@campus.haifa.ac.il;mhoush@univ.haifa.ac.il;elad.salomons@gmail.com
%% ####################################################################################################################
% This code requires:
% Developed under Matlab R2023a

%% Operational GHG emmision
OGHG = value(1.1 * sum(sum(sum(n_s_w_mat.*(E_GB+E_GP)))))*(1-(1+1^-6)^-100)/1^-6  % Operational GHG emmision is in (kg CO2-e) 

%% Present solution

% Define the decision variables and their values
decision_variables = {'Pipe diameter', 'Maximum tank volume', 'Head gain pump', 'Power of the pump', 'Flow rate of pump'}';
values = [double(d); double(V_max); double(H_p); double(P); double(Q_p)];
units = {'m', 'm^3', 'm', 'kW', 'm^3/hr'}';

% Concatenate the variables into a table
tab_w = table(decision_variables, values, units);
disp(tab_w)  % Display the table

% Define the decision variables and their values
decision_variables = {'Capacity of battery energy', ' Maximum battery discharging power','Minimum battery discharging power', 'Solar power Capacity'}';
values = [double(E_BST_max); double(P_BST_dmax); double(P_BST_cmax); double(max(P_S_std))];
units = {'kWh', 'kW', 'kW', 'kW'}';

% Concatenate the variables into a table
tab_p = table(decision_variables, values, units);

% Display the table
disp(tab_p)

%% Volume of the tank 

% Reshape matrix_1 to match the dimensions of matrix_2
V_tank_ini_1 = reshape(V_tank_ini, [1, size(V_tank_ini)]);

% Concatenate the matrices along the first dimension
concatenated_matrix = cat(1, V_tank_ini_1, V_tank);

vmin_column_1 = 300 * ones(25, 1); 
vmin_column_2 = 150 * ones(25, 1); 

resulting_matrix = [vmin_column_1, vmin_column_2, vmin_column_2, vmin_column_2]; 
vmax_column_1 =value(V_max) * ones(25, 1); 
resulting_matrix_1 = [vmax_column_1, vmax_column_1, vmax_column_1, vmax_column_1]; 

for w = 1:W
    figure('InvertHardcopy','off','Color',[1 1 1]);
    titles = {'(a) Summer','(b) Winter', '(c) Spring','(d) Autumn'};
    for s = 1:SEA
        subplot(2,2,s)
        title(titles(s),'FontWeight','bold', 'FontSize',8,'FontName','Times New Roman')
        hold on
        plot(0:24,resulting_matrix_1(:,s),'Color', 'k', linestyle = '-', LineWidth = 1)
        plot(0:24,resulting_matrix(:,s),'Color', 'k', linestyle = '--', LineWidth = 1)
        % Plot different color for specific ranges
        if s == 2
            if w == 1
                plot(0:17, value(concatenated_matrix(1:18,s,w)),'Color', [0 0.4470 0.7410], linestyle = '-', LineWidth = 1 ); 
                plot(17:22, value(concatenated_matrix(18:23,s,w)), 'Color','r',linestyle = '-', LineWidth = 1); 
                plot(22:24, value(concatenated_matrix(23:end,s,w)),'Color', [0 0.4470 0.7410], linestyle = '-',LineWidth = 1); 
                legend("V_{max}","V_{min}","Off peak_{WD}","On peak_{WD}","", "Orientation","horizontal", 'FontSize',7,...
                    'EdgeColor',[1 1 1],'FontWeight','bold', 'NumColumns', 2) 
            else
                plot(0:17, value(concatenated_matrix(1:18,s,w)), 'Color', [0 0.4470 0.7410], linestyle = '-', LineWidth = 1); 
                plot(17:22, value(concatenated_matrix(18:23,s,w)), 'r', linestyle = '-',LineWidth = 1); 
                plot(22:24, value(concatenated_matrix(23:end,s,w)), 'Color', [0 0.4470 0.7410],linestyle = '-', LineWidth = 1 ); 
                legend("V_{max}","V_{min}","Off peak_{WE}", "On peak_{WE}","Orientation","horizontal", 'FontSize',7,...
                    'EdgeColor',[1 1 1],'FontWeight','bold', 'NumColumns', 2) 
            end

        else
            if w == 1
                plot(0:17, value(concatenated_matrix(1:18,s,w)),'Color', [0 0.4470 0.7410], linestyle = '-',LineWidth = 1 ); 
                plot(17:23, value(concatenated_matrix(18:24,s,w)), 'r',linestyle = '-',LineWidth = 1); 
                plot(23:24, value(concatenated_matrix(24:end,s,w)),'Color', [0 0.4470 0.7410],  linestyle = '-',LineWidth = 1 );
            else
                plot(0:17, value(concatenated_matrix(1:18,s,w)),'Color', [0 0.4470 0.7410],  linestyle = '-',LineWidth = 1); 
                plot(17:23, value(concatenated_matrix(18:24,s,w)),'Color', [0 0.4470 0.7410], linestyle = '-',LineWidth = 1); 
                plot(23:24, value(concatenated_matrix(24:end,s,w)), 'Color', [0 0.4470 0.7410], linestyle = '-',LineWidth = 1); 
            end
        end
        if s==1|| s== 3
            ylabel(sprintf('Volume\n(m^3)'),'FontWeight','bold','FontName','Times New Roman','FontSize',9);
        end
        if s==3|| s==4
            xlabel('Time (hrs)','FontName','Times New Roman','FontSize',8,'FontWeight','bold');

        end
        ylim([0 2000])
        xlim([0 T])
        xticks(0:4:T)
        if s==1||s==2
            xticklabels([])
        end
        if s==2|| s== 4
            yticklabels([])
        end
        grid on
        box on
    end
end
 
 %% Pump Flow 

 for w= 1: W
     figure ('InvertHardcopy','off','Color',[1 1 1]);
     for s = 1:SEA
         subplot(2,2,s)
         title(titles(s),'FontWeight','bold', 'FontSize',8,'FontName','Times New Roman')
         xlim([0 T])
         xticks(0:4:T)
         grid on
         box on
         hold on
         stairs(flip([0 1:T]), flip([Q_dem(1,s,w); Q_dem(:,s,w)]),'LineStyle','-',LineWidth = 1);
         hold on
         stairs(flip([0 1:T]), flip([double(Q(1,s,w)); double(Q(:,s,w))]),'LineStyle','-',LineWidth = 1);
         if w == 1 && s == 4
             legend("Demand_{WD}", 'Pump flow_{WD}', "Orientation", "horizontal", 'FontSize', 7, ...
                 'EdgeColor', [1 1 1], 'FontWeight', 'bold', 'NumColumns', 2) 
         end
         if w == 2 && s == 4
             legend("Demand_{WE}", "Pump flow_{WE}", "Orientation", "horizontal", 'FontSize', 7, ...
                 'EdgeColor', [1 1 1], 'FontWeight', 'bold', 'NumColumns', 2) 
         end

         if s==1|| s== 3
             ylabel(sprintf('Flow \n(m^3 hr^{-1})'),'FontWeight','bold','FontName','Times New Roman','FontSize',9);
             ylim([0 600])
         end
         if s==2|| s== 4
             ylabel(sprintf(''),'FontWeight','bold','FontName','Times New Roman','FontSize',9);
             ylim([0 600])
             yticklabels([])
         end
         if s==3|| s==4
             xlabel('Time (hrs)','FontName','Times New Roman','FontSize',8,'FontWeight','bold');
         end
         if s==1||s==2
             xticklabels([])
         end
         hold off
     end
 end
 

%%  SOC 

% Reshape matrix_1 to match the dimensions of matrix_2
SOC_ini_1 = reshape(SOC_ini, [1, size(SOC_ini)]);

% Concatenate the matrices along the first dimension
concatenated_matrix_2 = cat(1, SOC_ini_1, SOC);

for w = 1:W
    figure ('InvertHardcopy','off','Color',[1 1 1]);
    for s= 1:SEA
        subplot(2,2,s)
        title(titles(s),'FontWeight','bold', 'FontSize',8,'FontName','Times New Roman')
        hold on
        plot(0:T, value(concatenated_matrix_2(:,s,w)),'LineStyle','-',LineWidth = 1);
        ylim([0 1])

        grid on
        box on
        if s==1|| s== 3
            ylabel('SOC','FontWeight','bold','FontName','Times New Roman','FontSize',8);
        end
        if s==3|| s==4
            xlabel('Time (hrs)','FontName','Times New Roman','FontSize',8,'FontWeight','bold');
        end
        xlim([0 T])
        xticks(0:4:T)
        if s==1||s==2
            xticklabels([])
        end
        if s==2|| s== 4
            yticklabels([])
        end
    end
    if w == 1 && s ==4
        legend("SOC_{WD}","Orientation","horizontal", 'FontSize',7,...
            'EdgeColor',[1 1 1],'FontWeight','bold') % Add a title and legend
    end
    if w == 2 && s == 4
        legend("SOC_{WE}","Orientation","horizontal", 'FontSize',7,...
            'EdgeColor',[1 1 1],'FontWeight','bold') % Add a title and legend
    end
    hold off
end

%% Energy to Pump 

for w = 1:W
    figure  ('InvertHardcopy','off','Color',[1 1 1]);
    for s = 1:SEA
        subplot(2,2,s);
        title(titles(s),'FontWeight','bold', 'FontSize',8,'FontName','Times New Roman')
        hold on;
        bar(1:24, [value(E_SP(:,s,w)*eta_inv), value(E_GP(:,s,w)), value(E_BP(:,s,w)*eta_inv*eta_d)]','stacked', 'BarWidth',1,'EdgeColor','none')
        grid on;
        box on;
        ylim([0 400])
        if s==1|| s== 3
            ylabel(sprintf('Energy to Pump \n(kWh)'),'FontWeight','bold','FontName','Times New Roman','FontSize',8);
        end
        if s==3|| s==4
            xlabel('Time (hrs)','FontName','Times New Roman','FontSize',8,'FontWeight','bold');
        end
        xlim([0 25])
        xticks(0:4:25)
        if s==1||s==2
            xticklabels([])
        end
        if s==2|| s== 4
            yticklabels([])
        end
    end
    legend("E_{SP WD}", "E_{GP WD}", "E_{BP WD}","Orientation","horizontal", 'FontSize',7,...
        'EdgeColor',[1 1 1],'FontWeight','bold')
    if w ==1
        legend("E_{SP WD}", "E_{GP WD}", "E_{BP WD}","Orientation","horizontal", 'FontSize',7,...
            'EdgeColor',[1 1 1],'FontWeight','bold')
    else
        legend("E_{SP WE}", "E_{GP WE}", "E_{BP WE}","Orientation","horizontal", 'FontSize',7,...
            'EdgeColor',[1 1 1],'FontWeight','bold')
    end

    hold off
end

%% Energy to grid

for w = 1:W
    figure('InvertHardcopy','off','Color',[1 1 1]);
    for s= 1:SEA
        subplot(2,2,s)
        title(titles(s),'FontWeight','bold', 'FontSize',8,'FontName','Times New Roman')
        hold on
        bar(1:T, [value(E_SG(:,s,w)*eta_inv)]','stacked','BarWidth',1,'EdgeColor','none')
        bar(1:T, [value(E_BG(:,s,w)*eta_inv*eta_d)]','stacked','BarWidth',1,'EdgeColor','none', FaceColor = [0.9290 0.6940 0.1250])
        grid on
        box on
        ylim([0 400])
        if s==1|| s== 3
            ylabel(sprintf("Selling to Grid \n(kWh)"),'FontWeight','bold','FontName','Times New Roman','FontSize',8);
        end
        if s==3|| s==4
            xlabel('Time (hrs)','FontName','Times New Roman','FontSize',8,'FontWeight','bold');
        end
        xlim([0 25])
        xticks(0:4:25)
        if s==1||s==2
            xticklabels([])
        end
        if s==2|| s== 4
            yticklabels([])
        end
    end
    if w ==1
        legend("E_{SG WD}","E_{BG WD}","Orientation","horizontal", 'FontSize',7,...
            'EdgeColor',[1 1 1],'FontWeight','bold')
    else
        legend("E_{SG WE}","E_{BG WE}","Orientation","horizontal", 'FontSize',7,...
            'EdgeColor',[1 1 1],'FontWeight','bold')
    end

    hold off
end
%% Energy to battery

for w = 1:W
    figure('InvertHardcopy','off','Color',[1 1 1]);
    for s= 1:SEA
        subplot(2,2,s)
        title(titles(s),'FontWeight','bold', 'FontSize',8,'FontName','Times New Roman')
        hold on
        bar(1:24, [value(E_SB(:,s,w)*eta_inv), value(E_GB(:,s,w)*eta_inv)]','stacked','BarWidth',1,'EdgeColor','none')
        grid on
        box on
        ylim([0 400])
        if s==1|| s== 3
            ylabel("to battery (kWh)",'FontWeight','bold','FontName','Times New Roman','FontSize',8);
        end
        if s==3|| s==4
            xlabel('Time (hrs)','FontName','Times New Roman','FontSize',8,'FontWeight','bold');
        end
        xlim([0 25])
        xticks(0:4:25)
        if s==1||s==2
            xticklabels([])
        end
        if s==2|| s== 4
            yticklabels([])
        end
    end
    if w == 1
        legend("E_{SB WD}","E_{GB WD}","Orientation","horizontal", 'FontSize',7,...
            'EdgeColor',[1 1 1],'FontWeight','bold')
    else
        legend("E_{SB WE}","E_{GB WE}","Orientation","horizontal", 'FontSize',7,...
            'EdgeColor',[1 1 1],'FontWeight','bold')

    end

    hold off
end