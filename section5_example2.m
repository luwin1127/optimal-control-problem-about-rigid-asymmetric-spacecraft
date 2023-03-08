%--------------------------------------------------------------------------
% ˵����Direct solution of nonlinear optimal control problems using 
%       quasilinearization and Chebyshev polynomials
%       pp.491-497 example 2
% ���ߣ�Lingwei Li
% ʱ�䣺ver 1.0 2023/03/6
%--------------------------------------------------------------------------
function section5_example2
%% �������
clc;clear;close all;

%% ���ó�ʼ����
% ״̬
w10 = .01;
w20 = .005;
w30 = .001;

w1_tf = 0;
w2_tf = 0;
w3_tf = 0;

% ����
I1 = 86.24;
I2 = 85.07;
I3 = 113.59;

%% α�׷����
nP = 75;                                % ���õ����
nS = 3;                                 % ״̬������
nU = 3;                                 % ����������

% α�׷����
tau = LGL_nodes(nP-1);                  % ��ʱ��������ɢ��[-1,1]��
D   = LGL_Dmatrix(tau);                 % ΢�־���
weights = LGL_weights(tau);             % ���õ�Ȩ��

%% �����㷨�ĳ�ֵ�²�
t0 = 0;                                 % ��ʼʱ��
T  = 100;                               % ��ֹʱ��
t  = tau2t(T,t0,tau);                   % �����ʱ������

% ��ǰ�²⼸��������
% �췽������
u10 = -.0087*ones(nP,1);
u20 = -.0043*ones(nP,1);
u30 = -.015*ones(nP,1);

% ������������
% ״̬��
w1_index = (1:nP)';
w2_index = (nP+1:2*nP)';
w3_index = (2*nP+1:3*nP)';
% ������
u1_index = (3*nP+1:4*nP)';
u2_index = (4*nP+1:5*nP)';
u3_index = (5*nP+1:6*nP)';

%% Ϊ R:max ���㷨���ó�ʼ����
x0 = zeros((nS+nU)*nP,1);
x0([w1_index,w2_index,w3_index]) = [w10*ones(nP,1);w20*ones(nP,1);w30*ones(nP,1)];
x0([w1_index(end),w2_index(end),w3_index(end)]) = [w1_tf;w2_tf;w3_tf];
x0([u1_index,u2_index,u3_index]) = [u10;u20;u30]; 

%% �㷨��ʼ
% ѭ������
options = optimset('Display','iter',...             % off/iter/notify-detailed
                   'Algorithm','sqp',...            % �㷨ѡ��
                   'TolCon',1e-6,...                % Լ��Υ���ȶ��ݲ����������Ĭ��ֵΪ 1e-6
                   'TolX',1e-6,...                  % ���������� x ����ֹ�ݲSQP�㷨��Ĭ��ֵΪ 1e-6
                   'TolFun',1e-6);               	% һ�������Ե���ֹ�ݲ����������Ĭ��ֵΪ 1e-6                  

tic;
%--- �㷨��ʼ ---%
res = fmincon(@objective,x0,[],[],[],[],[],[],@constraints,options);
%--- �㷨���� ---%
toc;

w1 = res(w1_index);
w2 = res(w2_index);
w3 = res(w3_index);
u1 = res(u1_index);
u2 = res(u2_index);
u3 = res(u3_index);
%% ��ͼ
% ��ʼ��ͼ
window_width = 500;
window_height = 416;

% ״̬��
k = 0;
figure('color',[1 1 1],'position',[300+k*window_width,300,window_width,window_height]);
plot(t, w1, '-',  'LineWidth',1.5);hold on;
plot(t, w2, '--', 'LineWidth',1.5);
plot(t, w3, ':',  'LineWidth',1.5);
legend('$\omega_{1}(t)$','$\omega_{2}(t)$','$\omega_{3}(t$)',...
       'Location','Best',...
       'Interpreter','Latex',...
       'FontSize',13);
xlabel('Time');
ylabel('State');
set(gca,'FontSize',13,'FontName','Times New Roman','LineWidth',1.5);

% ������
k = 1;
figure('color',[1 1 1],'position',[300+k*window_width,300,window_width,window_height]);
plot(t, u1, '-',  'LineWidth',1.5);hold on;
plot(t, u2, '--', 'LineWidth',1.5);
plot(t, u3, ':',  'LineWidth',1.5);
legend('$u_1(t)$','$u_2(t)$','$u_3(t)$',...
       'Location','Best',...
       'Interpreter','Latex',...
       'FontSize',13);
xlabel('Time');
ylabel('Control');
set(gca,'FontSize',13,'FontName','Times New Roman','LineWidth',1.5);

%% �Ӻ�������
function dJ = objective(x)
    % ��������ֵ
    u1 = x(u1_index);
    u2 = x(u2_index);
    u3 = x(u3_index);
    
    L = u1.^2 + u2.^2 + u3.^2;
    
    % Ŀ�꺯��
    dJ = (T-t0)/2 * dot(weights,L)/2;
end

function [c,ceq] = constraints(x)
    % ---------- ��ֵ ---------- %
    % ״̬����ֵ
    w1_con = x(w1_index);
    w2_con = x(w2_index);
    w3_con = x(w3_index);
    
    % ��������ֵ
    u1_con = x(u1_index);
    u2_con = x(u2_index);
    u3_con = x(u3_index);
    % ---------- ���� ---------- %
    
    % -------- ��ʽԼ�� -------- %
    ceq1 = w1_con(1) - w10;
    ceq2 = w2_con(1) - w20;
    ceq3 = w3_con(1) - w30;
    ceq4 = w1_con(end) - w1_tf;
    ceq5 = w2_con(end) - w2_tf;
    ceq6 = w3_con(end) - w3_tf;
    
    % ״̬����
    dw1_con = -(I3-I2)/I1 .* w2_con .* w3_con + u1_con ./ I1;
    dw2_con = -(I1-I3)/I2 .* w1_con .* w3_con + u2_con ./ I2;
    dw3_con = -(I2-I1)/I3 .* w1_con .* w2_con + u3_con ./ I3;
    
    Y = [w1_con,w2_con,w3_con];
    F = (T-t0)/2*real([dw1_con,dw2_con,dw3_con]);
    ceq7 = D*Y - F;
    ceq = [ceq1;ceq2;ceq3;ceq4;ceq5;ceq6;ceq7(:)];
    % ---------- ���� ---------- %
    
    % ------- ����ʽԼ�� ------- %
%     c = []; % case 1: no inequality constraints
    c = w1_con - (5.*1e-6.*t.^2 - 5.*1e-4.*t + .016); % case 2: inequality constraints on w1(t)
    % ---------- ���� ---------- %
end

function t = tau2t(Tf,t0,tau)
% �� tau ת��Ϊ t
% ���ա�����Եش����������滮���Ƶ�������2013 p79 ʽ(4.25)��д���´���
t = (Tf-t0)*tau/2+(Tf+t0)/2;
end
end