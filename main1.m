clc,clear
close all
%%
P_base = [0,0,100000];  %基准点位置
t_base = 0;
P_real = [-217960,112461,145800;
    266251,257393,111000;
    8056,154184,113100;
    -139758,181531,109300;
    177871,86010,105800;
     0,0,125000;
    ];

%P_real = [-217960,112461,145800;266251,257393,111000];
%%
mont = 1;   %实验次数
N_iter = 100;  % 迭代次数
c = 3e8;
N_sat = length(P_real); %卫星个数
N_sat = 6;
A_T = 1;   %时间误差方差ns
A_L = 10;  %位置误差方差cm
rmse_e_t = zeros(length(A_T),1);
error = zeros(mont,4*N_sat);
noise = 0e-2*randn(1,N_sat^2-N_sat);    %观测误差

% 预分配存储变量
T_real_sum = zeros(N_sat, 1);    % 用于累加 T_real
sT_sum = zeros(N_iter, N_sat);      % 用于累加 sT
sP_sum = zeros(N_iter, 3 * N_sat);  % 用于累加 sP

% 如果需要存储每次 Monte Carlo 仿真结果
T_real_store = zeros(mont, N_sat);         % 每次仿真的真实钟差
sT_store = zeros(mont, N_iter, N_sat);        % 每次仿真的 sT
sP_store = zeros(mont, N_iter, 3 * N_sat);    % 每次仿真的 sP


%%蒙特卡洛实验
for k = 1:mont
        %noise = 0e-2*randn(1,N_sat^2-N_sat);    %观测误差
        T_real = A_T*1e-9*c*randn(N_sat,1);     %每个卫星的钟差
        %T_real = 1e-9*c*[-0.3;-0.2;-0.5;0.1;0.2;0.3];
        L = A_L*1e-2*randn(N_sat,3);    %每个卫星的位置误差
        %L = A_L*1e-2*[0.1,0.1,0.1;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0];
        [T_find,P_find,sf,sT,sP] = ff1(P_real,P_base,T_real,L,noise,N_sat,N_iter);

        %%存储sf，sT和sK
        T_real_sum = T_real_sum + T_real;       % 累加真实钟差
        sT_sum = sT_sum + sT;                   % 累加 sT
        sP_sum = sP_sum + sP;                   % 累加 sP
        % 存储每次仿真结果（可选）
        T_real_store(k, :) = T_real';           % 存储真实钟差
        sT_store(k, :, :) = sT;                 % 存储 sT
        sP_store(k, :, :) = sP;                 % 存储 sP

end

T_real_avg = T_real_sum / mont;            % T_real 的 Monte Carlo 平均值
sT_k = sT_sum / mont;                      % sT 的 Monte Carlo 平均值
sP_k = sP_sum / mont;                      % sP 的 Monte Carlo 平均值

error_T = sT_k-repmat(T_real',N_iter,1); 
error_P = sP_k-reshape(P_real',1,[]);




%%画图
x = 1:N_iter;  % 横坐标：1 到 100
T_real_plot = repmat(T_real_avg, 1, N_iter);

figure
hold on
colors = lines(N_sat); % 获取颜色方案
for sat = 1:N_sat
    % 绘制 T_real_avg 的水平参考线
    plot(x, T_real_plot(sat, :), 'LineWidth', 0.5, 'Color', colors(sat, :),'DisplayName', ['T_{real,avg} - Sat ', num2str(sat)]);
    
    % 绘制 sT_k 的散点图
    scatter(x, sT_k(:, sat), 5, colors(sat, :),'DisplayName', ['sT_k - Sat ', num2str(sat)]);
end

%%绘制第一个卫星的坐标迭代散点图
figure
P_real_1 = P_real(3,:);
sP_k_1 = sP_k(:,7:9);
scatter3(P_real_1(1),P_real_1(2),P_real_1(3),"red");
hold on
scatter3(sP_k_1(:,1),sP_k_1(:,2),sP_k_1(:,3),"blue");


figure
plot(x,sum(error_P.^2,2))

figure
plot(x,sum(error_T.^2,2))










