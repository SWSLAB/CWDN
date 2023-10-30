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

% Formulate and solve the CWDN system.	

%%  Water decision variables

%d   -  Diameter of pipe (m).
%Q_p -  Flow of the pump (m^3/hr).
%H_p -  Head added by the pump (m).
%R   -  Resistance of the pipe.
%P   -  Power of the pump (kW).
%V_max -  Maximum volume of the tank (m^3).
%E_P   -  Energy of the pump (kWh).

V_tank_ini = sdpvar(SEA,W);   % Initial Volume of the tank (m^3). 
I = sdpvar(T,SEA,W);          % The time duration when the pump is ON.
V_tank = sdpvar(T,SEA,W);     % Volume of the tank (m^3).

%% Discretization options
Dopts= (SEA:W:30)*0.0254;
Popts=50:25:600;
Vopts=1000:500:10000;
% Build descritisation matrix for Q and H
if ~isfile('QQ_HH.mat')
    f1=@(Q_p,H_p,P)Q_p.*H_p*gama/(3600*eta)/1000-P;
    f2=@(Q_p,H_p,d)(2.7768*10^-6 *l)/(((C)^1.852)*(d^4.87))*Q_p.^1.852 + H_tank - H_res-H_p;
    for i=1:length(Popts)
        for j=1:length(Dopts)
            F=@(x)[f1(x(1),x(2),Popts(i));f2(x(1),x(2),Dopts(j))];
            xsol=fsolve(F,[0 0],optimoptions('fsolve','Display','none'));
            QQ(i,j)=xsol(1);
            HH(i,j)=xsol(2);
            % for plotting
            Qpvec=50:1000;
            Hpvec1=Popts(i)*1000*(3600*eta)/gama./Qpvec;
            Hpvec2=(2.7768*10^-6 *l)/(((C)^1.852)*(Dopts(j)^4.87))*Qpvec.^1.852 + H_tank - H_res;
            plot(Qpvec,Hpvec1,'-b')
            hold on
            plot(Qpvec,Hpvec2,'-r')
        end
    end
    ylim([0 3000])
    if any(any(abs(imag(QQ+HH))>1e-6))
        error('Complex Numbers')
    end
    save('QQ_HH.mat','QQ','HH')
else
    load ('QQ_HH.mat');
end

b_pipe=binvar(1,length(Dopts));
b_pump=binvar(1,length(Popts));
b_tank=binvar(1,length(Vopts));
Z=sdpvar(length(Popts),length(Dopts),'full');
for t=1:T
    for s=1:SEA
        for w=1:W
            Y{t,s,w}=sdpvar(length(Popts),length(Dopts),'full');
        end
    end
end

%% Water constraints

Q_p=sum(sum(Z.*QQ));
H_p=sum(sum(Z.*HH));
P=sum(b_pump.*Popts);
d=sum(b_pipe.*Dopts);
V_max=sum(b_tank.*Vopts);

C_aux1= [        0<=Z<=repmat(b_pipe,length(b_pump),1)];
C_aux1= [C_aux1; 0<=Z<=repmat(b_pump',1,length(b_pipe))];
C_aux1= [C_aux1; repmat(b_pump',1,length(b_pipe))+repmat(b_pipe,length(b_pump),1)-1<=Z<=1];

E_P=sdpvar(T,SEA,W); %just initilization to speed up the loop.
Q=sdpvar(T,SEA,W);   %just initilization to speed up the loop.
C_aux2=[];

for t=1:T
    for s=1:SEA
        for w=1:W
            E_P(t,s,w)=sum(sum(Y{t,s,w}.*QQ.*HH*gama/(3600*eta)/1000)*dt);
            Q(t,s,w)=sum(sum(Y{t,s,w}.*QQ));
            C_aux2=[C_aux2; 0<=Y{t,s,w}<=I(t,s,w)];
            C_aux2=[C_aux2; I(t,s,w)-dt*(1-Z)<=Y{t,s,w}<=dt*Z];
        end
    end
end

for s = 1:SEA 
    for w = 1:W
        V_tank(:,s,w) = V_tank_ini(s,w) + cumsum(Q(:,s,w) *dt) - cumsum(Q_dem(:,s,w)*dt);
    end 
end

C1_W =[sum(b_pump)==1;sum(b_pipe)==1;sum(b_tank)==1];
C2_W = [0 <= I <= dt];
C3_W = [];
C4_W = [];
for s = 1:SEA
    for w = 1:W
        C3_W = [C3_W,  V_min(s,w) <= V_tank(:,s,w) <= V_max];
        C4_W = [C4_W, V_tank_ini(s,w) == V_tank(T,s,w)] ;
    end
end

Const_W = [C1_W, C2_W, C3_W, C4_W,C_aux1,C_aux2];

%% Water objective function

CPC = @(d)((32.598+0.11*(d*1000)+0.00053*(d*1000).^2)*l)*1.1*Crf_pipe;                  % Cost of the Pipe ($).
CPUC = @(P)(11603*(P.^0.47) + 42853*(P.^0.41))*1.1*Crf_pump;                            % Cost of the Pump ($).
CTC = @(V_max)(11709*(V_max.^0.3) + 2388*(V_max.^0.62))*1.1*Crf_tank;                   % Cost of the Tank ($).

obj_W = sum(CPUC(Popts).*b_pump) + sum(CTC(Vopts).*b_tank) + sum(CPC(Dopts).*b_pipe);   % Total cost of WDN.

%%  Power decision variables

%E_S        -   Energy of the Solar (kWh).
%E_GP       -   Energy from Grid to pump (kWh).
%E_SG       -   Energy from Solar to Grid (kWh).
%P_S_std    -   Power of the solar plant (kW).
%E_BST_max  -   Energy of the battery required (kWh).
%P_BST_dmax -   Power of the battery discharging(kW).
%P_BST_cmax -   Power of the battery charging(kW).

N_pan = intvar(1, 1);              % No. of panels in the solar plant.
E_BST_ini= sdpvar(SEA,W);          % Initial energy in the battery (kWh).
E_BST = sdpvar(T,SEA,W);           % Energy of the battery (kWh).
E_SP = sdpvar(T,SEA,W);            % Energy from Solar to pump (kWh),
E_BP = sdpvar(T,SEA,W);            % Energy from Battery to pump (kWh).
E_BG = sdpvar(T,SEA,W);            % Energy from Battery to Grid (kWh).
E_SB = sdpvar(T,SEA,W);            % Energy from Solar to Battery (kWh).
E_GB = sdpvar(T,SEA,W);            % Energy from Grid to Battery (kWh).
SOC = sdpvar(T,SEA,W);             % State of charge of battery.

%% Discretization options
E_BST_max_opt=0:50:1000;
b_bat=binvar(1,length(E_BST_max_opt));
E_BST_max=sum(b_bat.*E_BST_max_opt);
  
%% Power constraints

P_S_std = A_p*eta_pan *IR_std*DF* N_pan; 
E_S = (P_S_std.*(IR_t/IR_std)).* dt;
E_SG = E_S - E_SP - E_SB;
E_GP = E_P - E_SP*eta_inv - (E_BP*eta_inv)*eta_d;
P_BST_cmax = E_BST_max/C_c_time;
P_BST_dmax = E_BST_max/C_d_time;

for s = 1:SEA
    for w = 1:W
        E_BST(:, s, w) = E_BST_ini(s, w) + cumsum((E_GB(:, s, w) * eta_inv + E_SB(:, s, w)) * eta_c) - cumsum(E_BG(:, s, w) + E_BP(:, s, w));
        SOC(:,s,w) = E_BST(:,s,w)/E_BST_max;
        
    end
end

SOC_ini = E_BST_ini/ E_BST_max;

C1_P = (E_BP + E_BG)*eta_d <= P_BST_dmax*dt;
C2_P = (E_GB*eta_inv + E_SB)* eta_c <= P_BST_cmax*dt;
C3_P = [A_p*N_pan <= Amax];

C4_P = [];
C5_P = [];
for s = 1:SEA
    for w = 1:W
        C4_P = [C4_P, SOC_min*E_BST_max <= E_BST(:,s,w) <= SOC_max*E_BST_max];
        C5_P  = [C5_P,  E_BST_ini(s, w) == E_BST(T, s, w)];
    end
end

C6_P = [E_S >= 0; E_SG >= 0 ;E_SB >= 0;E_SP >= 0; E_BG >= 0;E_BP >= 0;N_pan>=0;E_BST_max>=0;P_BST_cmax>=0;P_BST_dmax>=0;E_BST>=0; E_GB>=0;E_GP >=0; 0 <= E_BST_ini <= E_BST_max];

Const_P = [C1_P, C2_P, C3_P, C4_P, C5_P, C6_P];

%% Power objective function

COS= @(A,P) Crf_solar*(382/3.5*A) + 0.0165*P;                           % Cost of the Solar ($).
COB = @(E_BST_max) Crf_bat*(150*(E_BST_max))+ 0.0165*(P_BST_dmax);      % Cost of the Battery source ($).

obj_P = COS(A_p*N_pan,P_S_std) + COB(E_BST_max) - sum(sum(sum((E_S_tariff.*(E_SG + (E_BG*eta_d))*eta_inv).*n_s_w_mat))) + sum(sum(sum(E_tariff.*(E_GB + E_GP).*n_s_w_mat))) ; % Total Cost of PDN.

%% total Objective function
obj_eps=(sum(E_BST_ini(:))+ sum(V_tank_ini(:))+sum((E_GP +E_SP + E_BP + E_GB + E_SB + E_BG + E_SG).* repmat((1:24)'*1/24,1,4,2),'all'))
obj = obj_W + obj_P+1e-5*obj_eps;

%% Define and solve the optimization problem

ops = sdpsettings('verbose', 5,'solver','gurobi');
ops.gurobi.MIPGap=1e-6;
sol = optimize([Const_P,Const_W], obj, ops);

% Extract and display the optimal solution:
if sol.problem == 0
    fprintf('Optimal solution found!\n');
    cost_for_Power_and_water_distribution_network = value(obj)
    Cost_of_WDN = value(obj_W)
    Cost_of_PDN = value(obj_P)
else
    error('Problem solving optimization problem!');
end




  