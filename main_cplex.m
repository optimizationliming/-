  %% ����������ϵͳ��΢���������Ż�_����ï
% ժҪ���ڵ�ǰ��Դ������Ѹ�ٷ�չ��������ϵ�ս����ܵĻ����£�������ڵ������ϵ��ȵ���������΢���������Ż�
% ģ�͡��ۺ����ڴ������ԡ���ʱ��ۡ����ȸ�����ֲ�ʽ��Դ��ʱ���������԰�������������ء��ȵ�����ϵͳ����
% ��¯��ȼ�ϵ�غʹ���ϵͳ�Ĳ�����΢����Ϊ��������Cplex�Ż������õ��������ڸ�΢��Դ��ѳ����������гɱ���
% �������ֳ������ȵ��ȷ�ʽ���бȽϡ������������������ϵ���ģ����ʵ�ֵ���ͳһЭ�����Ȳ�����΢�������гɱ���
% ��ģ�Ϳ�Ϊ����֮����Դ�������滮��Ӫ�ṩ�ο���

clc; 
clear; 
close all; 
warning off;

%% ������߱���
% �������
parameter_cplex;
Ppv_1 = sdpvar(1,24);
Pwt_1 = sdpvar(1,24);
% ���豸
Pmt = sdpvar(1,24);            % ȼ���ֻ��繦��
Peb = sdpvar(1,24);            % ���¯�繦�� 
Pfc = sdpvar(1,24);              % ȼ�ϵ�ص繦��
Pgrid = sdpvar(1,24);             % ��������                   >0 ��磬<0 ����
Eb= sdpvar(1,24);               % �索�ܵ�����
Pbch = sdpvar(1,24);          % �索�ܳ�繦��
Pbdis = sdpvar(1,24);         % �索�ܷŵ繦��
% ���豸
Pmth = sdpvar(1,24);          % ȼ���ֻ��ȹ���
Pheb = sdpvar(1,24);          % ���¯�ȹ���
Eh= sdpvar(1,24);               % �ȴ��ܵ�����
Phch = sdpvar(1,24);          % �ȴ��ܴ��ȹ���
Phdis = sdpvar(1,24);         % �ȴ��ܷ��ȹ���
% ��������
Ubch = binvar(1,24);   % ��س��״̬��1��ʾ���
Ubdis = binvar(1,24);   % ��طŵ�״̬��1��ʾ�ŵ�
Uhch = binvar(1,24);   % �ȴ��ܴ���״̬��1��ʾ����
Uhdis = binvar(1,24);   % �ȴ��ܷ���״̬��1��ʾ����
% ��ʼ��
objective=0;
constraint=[];

%% д���������Լ��
constraint=[constraint,0.8*Ppv<=Ppv_1<=1*Ppv];
constraint=[constraint,0.8*Pwt<=Pwt_1<=1*Pwt];

%% ȼ���ֻ�ģ��
constraint=[constraint, Pmth == ((Pmt.*(1 - eta_mt - eta_l))./eta_mt) * eta_h * Coph];
%% ���¯ģ��
constraint=[constraint, Pheb == Peb * eta_ah];
%% �繦��ƽ��Լ��
constraint = [constraint, Pmt + Pfc + Pgrid + Pwt_1 + Ppv_1+ Pbdis == Pl + Pbch + Peb];
%% �ȹ���ƽ��Լ��
constraint = [constraint, Pmth + Pheb + Phdis == Phe + Phch];
%% �索��Լ��
% 1. �索�ܲ���ʽԼ��
constraint=[constraint, Ubch*Pbch_min <= Pbch <= Ubch*Pbch_max];
constraint=[constraint, Ubdis*Pbdis_min <= Pbdis <= Ubdis*Pbdis_max];
constraint=[constraint, Ubch + Ubdis <= 1];
constraint=[constraint, Eb_min <= Eb <= Eb_max];
% 2. �索�ܵ�ʽԼ��
for t=1:24
     if t==1
        constraint=[constraint, Eb(t) == Eb_init*(1-tau_b) + Pbch(t)*eta_bch - Pbdis(t)/eta_bdis];
    else
        constraint=[constraint, Eb(t) == Eb(t-1)*(1-tau_b) + Pbch(t)*eta_bch - Pbdis(t)/eta_bdis];
    end
end
% 3. �索��ʼĩ���Լ��
constraint = [constraint, Eb_init == Eb(24)];
%% �ȴ���Լ��
% 1. �ȴ��ܲ���ʽԼ��
constraint=[constraint, Uhch*Phch_min <= Phch <= Uhch*Phch_max];
constraint=[constraint, Uhdis*Phdis_min <= Phdis <= Uhdis*Phdis_max];
constraint=[constraint, Uhch + Uhdis <= 1];
constraint=[constraint, Eh_min <= Eh <= Eh_max];
% 2. �ȴ��ܵ�ʽԼ��
for t=1:24
     if t==1
        constraint=[constraint, Eh(t) == Eh_init*(1-tau_h) + Phch(t)*eta_hch - Phdis(t)/eta_hdis];
    else
        constraint=[constraint, Eh(t) == Eh(t-1)*(1-tau_h) + Phch(t)*eta_hch - Phdis(t)/eta_hdis];
    end
end
% 3. �ȴ���ʼĩ���Լ��
constraint = [constraint, Eh_init == Eh(24)];
 %% ���߱����߽�Լ��
constraint=[constraint, Pmt_min <= Pmt <= Pmt_max];
constraint=[constraint, Peb_min <= Peb <= Peb_max];
constraint=[constraint, Pfc_min <= Pfc <= Pfc_max];
constraint=[constraint, Pgrid_min <= Pgrid <= Pgrid_max];
%% �����豸����Լ��
for t=2:24
    constraint=[constraint, Rmt_down <= Pmt(t) - Pmt(t-1) <= Rmt_up];
    constraint=[constraint, Rfc_down <= Pfc(t) - Pfc(t-1) <= Rfc_up];
    constraint=[constraint, Reb_down <= Peb(t) - Peb(t-1) <= Reb_up];
end

%% Ŀ�꺯�����յ��ȳɱ���С
for t=1:24
    objective = objective + Cgas*( (Pmt(t)/(L_gas*eta_mt(t))) + (Pfc(t)/(L_gas*eta_fc(t))) ) ...  % ȼ�ϳɱ�
                         + 1/2*(Crb(t)+Crs(t))*Pgrid(t) + 1/2*(Crb(t)-Crs(t))*abs(Pgrid(t)) ...                  % �����ɱ�
                         + Cm_mt*Pmt(t) +Cm_fc*Pfc(t)  + Cm_eb*Peb(t) + Cm_wt*Ppv_1(t)...          % ά���ɱ�
                         + Cm_pv*Pwt_1(t) + Cm_Eb*(Pbch(t)+Pbdis(t)) + Cm_Eh*(Phch(t)+Phdis(t))...                     
                          - Che*Phe(t);                                                                                                         % ��������
end

%% ���
option = sdpsettings('verbose',1,'solver','cplex');
Diagnostic_information = optimize(constraint,objective,option)
if Diagnostic_information.problem == 0
    disp('Diagnostic_information����ģ���н⣬���Ž�����');
    Cost_op = value(objective)
%     check(constraint)
else
    disp('Diagnostic_information����ģ���޽⣬����Լ���Ƿ��ͻ������Ƿ����');
    check(constraint)
end

%% �Ѹ��������ÿСʱ������yalmip��ȡ����
Pmt_in = value(Pmt);                                           % ȫΪ+
Pfc_in = value(Pfc);                                              % ȫΪ+
Peb_out = value(Peb);                                          % ȫΪ+
Pmth_in = value(Pmth);                                       % ȫΪ+
Pheb_in = value(Pheb);                                        % ȫΪ+
Pgrid_result = value(Pgrid);                                  % �����и�
Pb_result= value(Pbdis-Pbch);                              % �����и�  �索��
Ph_result = value(Phdis-Phch);  % �����и� �ȴ���
bSOC = value((Eb/150)*100);
% Pwt_1_in = value(Pwt_1);
% Ppv_1_in = value(Ppv_1);
hSOC = value((Eh/100)*100);

%% ��ƽ��ͼ
figure(1)
hold on 
% ����
Pgrid_in = max(Pgrid_result, 0);
Pb_in = max(Pb_result, 0);
electricity_in = [Pmt_in', Pfc_in', Pgrid_in', Pb_in'];
% �õ�
Pgrid_out = min(Pgrid_result, 0);
Pb_out = min(Pb_result, 0);
electricity_out = [Pgrid_out', Pb_out', (-Peb_out)'];
% ��ͼ
bar(electricity_in, 'stack');
bar(electricity_out, 'stack');
plot(Pl-Ppv-Pwt, 'ok-');
plot(bSOC, '.r:');
title('electric power balance');
legend('Pmt','Pfc', 'Pex buy', 'Pbdis', 'Pex sell', 'Pbch', 'Peb', '���縺��', 'SOC_��');
grid on
hold off

%% ��ƽ��ͼ
figure(2)
hold on 
% ����
Ph_in = max(Ph_result, 0);
Thermal_in = [Pmth_in', Pheb_in', Ph_in'];
% ����
Ph_out = min(Ph_result, 0);
Thermal_out = Ph_out';
% ��ͼ
bar(Thermal_in, 'stack');
bar(Thermal_out, 'stack');
plot(Phe, 'sk-');
plot(hSOC, '.r:');
title('thermal power balance');
legend('Pmth', 'Pheb', 'Phdis', 'Phch', '�ȸ���', 'SOC_��');
grid on
hold off
