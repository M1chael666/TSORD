function [T_find,P_find,sf,sT,sP] = ff1(P_real,P_base,T_real,L,noise, N_sat,N_iter)

%%各卫星到参考点的距离
t_N0_real = zeros(1,N_sat);
for m = 1:N_sat
    t_N0_real(1,m) = sqrt(sum((P_real(m,:)-P_base).^2));
end

%%各卫星之间距离
t_inter_real = zeros(N_sat,N_sat);
for m =1:N_sat
    for n =1:N_sat
        t_inter_real(m,n) = sqrt(sum((P_real(m,:)-P_real(n,:)).^2));
    end
end

%%观测值
t_inter_estimate = zeros(N_sat,N_sat);
i=0;
for m = 1:N_sat
    for n =1:N_sat
        if m==n
            t_inter_estimate(m,n)= t_inter_estimate(m,n)+0;
        else
            i=i+1;
            t_inter_estimate(m,n) = T_real(n)-(t_inter_real(m,n)-t_N0_real(1,m))+noise(i);
        end
    end
end

%%认为初始无钟差
T_0 = zeros(1,N_sat); %行向量
%位置初值
P_0 = P_real+L;


%%梯度迭代
for i=1:N_iter
    [dj,f] = partial1(T_0,P_0,P_base,t_inter_estimate,N_sat);
    f_obj_val = f;
    sf(i) = f;
    sT(i,:)= T_0;
    sP(i,:)= reshape(P_0', 1, []);

    %%回溯线搜索
    stepsize = 1;   %步长
    T_1 = T_0;
    P_1 = P_0;
    dj0 = dj;
    f_obj_val2 = f;
    beta = 0.2;  %衰减因子
    alpha = 1e-9;
    ttk = 0;

    while (f_obj_val-f_obj_val2)<alpha*stepsize*norm(dj).^2
        stepsize = stepsize*beta;
        T_0 = T_1-stepsize*dj(1:N_sat);
        P_0 = P_1-stepsize*reshape(dj(N_sat+1:end), 3, []).';

        [dj,f] = partial1(T_0,P_0,P_base,t_inter_estimate,N_sat);
        f_obj_val2 = f;
        ttk = ttk+1;
        if ttk>=150
            break;
        end
    end
    
    T_0 = T_1-stepsize*dj0(1:N_sat);
    P_0 = P_1-stepsize*reshape(dj0(N_sat+1:end), 3, []).';

end

T_find = T_real.';
P_find = P_0-P_real;










