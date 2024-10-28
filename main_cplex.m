  %% 含电热联合系统的微电网运行优化_李正茂
% 摘要：在当前能源互联网迅速发展及电热联系日渐紧密的环境下，提出基于电热联合调度的区域并网型微电网运行优化
% 模型。综合网内储能特性、分时电价、电热负荷与分布式电源的时序特征，以包含风机、光伏电池、热电联产系统、电
% 锅炉、燃料电池和储能系统的并网型微电网为例，采用Cplex优化软件求得调度周期内各微电源最佳出力及总运行成本，
% 并与两种常见电热调度方式进行比较。仿真算例表明：联合调度模型能实现电热统一协调调度并降低微电网运行成本。
% 该模型可为电热之间能源互联及规划运营提供参考。

clc; 
clear; 
close all; 
warning off;

%% 定义决策变量
% 输入参数
parameter_cplex;
Ppv_1 = sdpvar(1,24);
Pwt_1 = sdpvar(1,24);
% 电设备
Pmt = sdpvar(1,24);            % 燃气轮机电功率
Peb = sdpvar(1,24);            % 电锅炉电功率 
Pfc = sdpvar(1,24);              % 燃料电池电功率
Pgrid = sdpvar(1,24);             % 电网功率                   >0 买电，<0 卖电
Eb= sdpvar(1,24);               % 电储能的容量
Pbch = sdpvar(1,24);          % 电储能充电功率
Pbdis = sdpvar(1,24);         % 电储能放电功率
% 热设备
Pmth = sdpvar(1,24);          % 燃气轮机热功率
Pheb = sdpvar(1,24);          % 电锅炉热功率
Eh= sdpvar(1,24);               % 热储能的容量
Phch = sdpvar(1,24);          % 热储能储热功率
Phdis = sdpvar(1,24);         % 热储能放热功率
% 辅助变量
Ubch = binvar(1,24);   % 电池充电状态，1表示充电
Ubdis = binvar(1,24);   % 电池放电状态，1表示放电
Uhch = binvar(1,24);   % 热储能储热状态，1表示储热
Uhdis = binvar(1,24);   % 热储能放热状态，1表示放热
% 初始化
objective=0;
constraint=[];

%% 写个风光区间约束
constraint=[constraint,0.8*Ppv<=Ppv_1<=1*Ppv];
constraint=[constraint,0.8*Pwt<=Pwt_1<=1*Pwt];

%% 燃气轮机模型
constraint=[constraint, Pmth == ((Pmt.*(1 - eta_mt - eta_l))./eta_mt) * eta_h * Coph];
%% 电锅炉模型
constraint=[constraint, Pheb == Peb * eta_ah];
%% 电功率平衡约束
constraint = [constraint, Pmt + Pfc + Pgrid + Pwt_1 + Ppv_1+ Pbdis == Pl + Pbch + Peb];
%% 热功率平衡约束
constraint = [constraint, Pmth + Pheb + Phdis == Phe + Phch];
%% 电储能约束
% 1. 电储能不等式约束
constraint=[constraint, Ubch*Pbch_min <= Pbch <= Ubch*Pbch_max];
constraint=[constraint, Ubdis*Pbdis_min <= Pbdis <= Ubdis*Pbdis_max];
constraint=[constraint, Ubch + Ubdis <= 1];
constraint=[constraint, Eb_min <= Eb <= Eb_max];
% 2. 电储能等式约束
for t=1:24
     if t==1
        constraint=[constraint, Eb(t) == Eb_init*(1-tau_b) + Pbch(t)*eta_bch - Pbdis(t)/eta_bdis];
    else
        constraint=[constraint, Eb(t) == Eb(t-1)*(1-tau_b) + Pbch(t)*eta_bch - Pbdis(t)/eta_bdis];
    end
end
% 3. 电储能始末相等约束
constraint = [constraint, Eb_init == Eb(24)];
%% 热储能约束
% 1. 热储能不等式约束
constraint=[constraint, Uhch*Phch_min <= Phch <= Uhch*Phch_max];
constraint=[constraint, Uhdis*Phdis_min <= Phdis <= Uhdis*Phdis_max];
constraint=[constraint, Uhch + Uhdis <= 1];
constraint=[constraint, Eh_min <= Eh <= Eh_max];
% 2. 热储能等式约束
for t=1:24
     if t==1
        constraint=[constraint, Eh(t) == Eh_init*(1-tau_h) + Phch(t)*eta_hch - Phdis(t)/eta_hdis];
    else
        constraint=[constraint, Eh(t) == Eh(t-1)*(1-tau_h) + Phch(t)*eta_hch - Phdis(t)/eta_hdis];
    end
end
% 3. 热储能始末相等约束
constraint = [constraint, Eh_init == Eh(24)];
 %% 决策变量边界约束
constraint=[constraint, Pmt_min <= Pmt <= Pmt_max];
constraint=[constraint, Peb_min <= Peb <= Peb_max];
constraint=[constraint, Pfc_min <= Pfc <= Pfc_max];
constraint=[constraint, Pgrid_min <= Pgrid <= Pgrid_max];
%% 部分设备爬坡约束
for t=2:24
    constraint=[constraint, Rmt_down <= Pmt(t) - Pmt(t-1) <= Rmt_up];
    constraint=[constraint, Rfc_down <= Pfc(t) - Pfc(t-1) <= Rfc_up];
    constraint=[constraint, Reb_down <= Peb(t) - Peb(t-1) <= Reb_up];
end

%% 目标函数：日调度成本最小
for t=1:24
    objective = objective + Cgas*( (Pmt(t)/(L_gas*eta_mt(t))) + (Pfc(t)/(L_gas*eta_fc(t))) ) ...  % 燃料成本
                         + 1/2*(Crb(t)+Crs(t))*Pgrid(t) + 1/2*(Crb(t)-Crs(t))*abs(Pgrid(t)) ...                  % 电网成本
                         + Cm_mt*Pmt(t) +Cm_fc*Pfc(t)  + Cm_eb*Peb(t) + Cm_wt*Ppv_1(t)...          % 维护成本
                         + Cm_pv*Pwt_1(t) + Cm_Eb*(Pbch(t)+Pbdis(t)) + Cm_Eh*(Phch(t)+Phdis(t))...                     
                          - Che*Phe(t);                                                                                                         % 制热收益
end

%% 求解
option = sdpsettings('verbose',1,'solver','cplex');
Diagnostic_information = optimize(constraint,objective,option)
if Diagnostic_information.problem == 0
    disp('Diagnostic_information：该模型有解，最优解如下');
    Cost_op = value(objective)
%     check(constraint)
else
    disp('Diagnostic_information：该模型无解，请检查约束是否冲突或参数是否合理');
    check(constraint)
end

%% 把各个机组的每小时出力从yalmip中取出来
Pmt_in = value(Pmt);                                           % 全为+
Pfc_in = value(Pfc);                                              % 全为+
Peb_out = value(Peb);                                          % 全为+
Pmth_in = value(Pmth);                                       % 全为+
Pheb_in = value(Pheb);                                        % 全为+
Pgrid_result = value(Pgrid);                                  % 有正有负
Pb_result= value(Pbdis-Pbch);                              % 有正有负  电储能
Ph_result = value(Phdis-Phch);  % 有正有负 热储能
bSOC = value((Eb/150)*100);
% Pwt_1_in = value(Pwt_1);
% Ppv_1_in = value(Ppv_1);
hSOC = value((Eh/100)*100);

%% 电平衡图
figure(1)
hold on 
% 发电
Pgrid_in = max(Pgrid_result, 0);
Pb_in = max(Pb_result, 0);
electricity_in = [Pmt_in', Pfc_in', Pgrid_in', Pb_in'];
% 用电
Pgrid_out = min(Pgrid_result, 0);
Pb_out = min(Pb_result, 0);
electricity_out = [Pgrid_out', Pb_out', (-Peb_out)'];
% 画图
bar(electricity_in, 'stack');
bar(electricity_out, 'stack');
plot(Pl-Ppv-Pwt, 'ok-');
plot(bSOC, '.r:');
title('electric power balance');
legend('Pmt','Pfc', 'Pex buy', 'Pbdis', 'Pex sell', 'Pbch', 'Peb', '净电负荷', 'SOC_电');
grid on
hold off

%% 热平衡图
figure(2)
hold on 
% 制热
Ph_in = max(Ph_result, 0);
Thermal_in = [Pmth_in', Pheb_in', Ph_in'];
% 耗热
Ph_out = min(Ph_result, 0);
Thermal_out = Ph_out';
% 画图
bar(Thermal_in, 'stack');
bar(Thermal_out, 'stack');
plot(Phe, 'sk-');
plot(hSOC, '.r:');
title('thermal power balance');
legend('Pmth', 'Pheb', 'Phdis', 'Phch', '热负荷', 'SOC_热');
grid on
hold off
