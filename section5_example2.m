%--------------------------------------------------------------------------
% 说明：Direct solution of nonlinear optimal control problems using 
%       quasilinearization and Chebyshev polynomials
%       pp.491-497 example 2
% 作者：Lingwei Li
% 时间：ver 1.0 2023/03/6
%--------------------------------------------------------------------------
function section5_example2
%% 清除变量
clc;clear;close all;

%% 设置初始条件
% 状态
w10 = .01;
w20 = .005;
w30 = .001;

w1_tf = 0;
w2_tf = 0;
w3_tf = 0;

% 参数
I1 = 86.24;
I2 = 85.07;
I3 = 113.59;

%% 伪谱法配点
nP = 75;                                % 配置点个数
nS = 3;                                 % 状态量个数
nU = 3;                                 % 控制量个数

% 伪谱法配点
tau = LGL_nodes(nP-1);                  % 把时间区间离散到[-1,1]中
D   = LGL_Dmatrix(tau);                 % 微分矩阵
weights = LGL_weights(tau);             % 配置点权重

%% 处理算法的初值猜测
t0 = 0;                                 % 起始时间
T  = 100;                               % 终止时间
t  = tau2t(T,t0,tau);                   % 配点后的时间区段

% 提前猜测几个控制量
% 红方控制量
u10 = -.0087*ones(nP,1);
u20 = -.0043*ones(nP,1);
u30 = -.015*ones(nP,1);

% 建立变量索引
% 状态量
w1_index = (1:nP)';
w2_index = (nP+1:2*nP)';
w3_index = (2*nP+1:3*nP)';
% 控制量
u1_index = (3*nP+1:4*nP)';
u2_index = (4*nP+1:5*nP)';
u3_index = (5*nP+1:6*nP)';

%% 为 R:max 的算法设置初始条件
x0 = zeros((nS+nU)*nP,1);
x0([w1_index,w2_index,w3_index]) = [w10*ones(nP,1);w20*ones(nP,1);w30*ones(nP,1)];
x0([w1_index(end),w2_index(end),w3_index(end)]) = [w1_tf;w2_tf;w3_tf];
x0([u1_index,u2_index,u3_index]) = [u10;u20;u30]; 

%% 算法开始
% 循环参数
options = optimset('Display','iter',...             % off/iter/notify-detailed
                   'Algorithm','sqp',...            % 算法选择
                   'TolCon',1e-6,...                % 约束违反度度容差（正标量），默认值为 1e-6
                   'TolX',1e-6,...                  % 关于正标量 x 的终止容差，SQP算法的默认值为 1e-6
                   'TolFun',1e-6);               	% 一阶最优性的终止容差（正标量），默认值为 1e-6                  

tic;
%--- 算法开始 ---%
res = fmincon(@objective,x0,[],[],[],[],[],[],@constraints,options);
%--- 算法结束 ---%
toc;

w1 = res(w1_index);
w2 = res(w2_index);
w3 = res(w3_index);
u1 = res(u1_index);
u2 = res(u2_index);
u3 = res(u3_index);
%% 画图
% 开始画图
window_width = 500;
window_height = 416;

% 状态量
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

% 控制量
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

%% 子函数部分
function dJ = objective(x)
    % 控制量赋值
    u1 = x(u1_index);
    u2 = x(u2_index);
    u3 = x(u3_index);
    
    L = u1.^2 + u2.^2 + u3.^2;
    
    % 目标函数
    dJ = (T-t0)/2 * dot(weights,L)/2;
end

function [c,ceq] = constraints(x)
    % ---------- 赋值 ---------- %
    % 状态量赋值
    w1_con = x(w1_index);
    w2_con = x(w2_index);
    w3_con = x(w3_index);
    
    % 控制量赋值
    u1_con = x(u1_index);
    u2_con = x(u2_index);
    u3_con = x(u3_index);
    % ---------- 结束 ---------- %
    
    % -------- 等式约束 -------- %
    ceq1 = w1_con(1) - w10;
    ceq2 = w2_con(1) - w20;
    ceq3 = w3_con(1) - w30;
    ceq4 = w1_con(end) - w1_tf;
    ceq5 = w2_con(end) - w2_tf;
    ceq6 = w3_con(end) - w3_tf;
    
    % 状态方程
    dw1_con = -(I3-I2)/I1 .* w2_con .* w3_con + u1_con ./ I1;
    dw2_con = -(I1-I3)/I2 .* w1_con .* w3_con + u2_con ./ I2;
    dw3_con = -(I2-I1)/I3 .* w1_con .* w2_con + u3_con ./ I3;
    
    Y = [w1_con,w2_con,w3_con];
    F = (T-t0)/2*real([dw1_con,dw2_con,dw3_con]);
    ceq7 = D*Y - F;
    ceq = [ceq1;ceq2;ceq3;ceq4;ceq5;ceq6;ceq7(:)];
    % ---------- 结束 ---------- %
    
    % ------- 不等式约束 ------- %
%     c = []; % case 1: no inequality constraints
    c = w1_con - (5.*1e-6.*t.^2 - 5.*1e-4.*t + .016); % case 2: inequality constraints on w1(t)
    % ---------- 结束 ---------- %
end

function t = tau2t(Tf,t0,tau)
% 把 tau 转换为 t
% 按照《天基对地打击武器轨道规划与制导技术》2013 p79 式(4.25)编写如下代码
t = (Tf-t0)*tau/2+(Tf+t0)/2;
end
end