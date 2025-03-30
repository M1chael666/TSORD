function [dj,f] = partial1(T_0,P_0,P_base,t_inter_estimate,N_sat)

%%各卫星到参考点
t_N0 = zeros(1,N_sat);
for m = 1:N_sat
    t_N0(1,m) = sqrt(sum((P_0(m,:)-P_base).^2));
end

%%各卫星之间
t_inter = zeros(N_sat,N_sat);
for m = 1:N_sat
    for n = 1:N_sat
        t_inter(m,n) = sqrt(sum((P_0(m,:)-P_0(n,:)).^2));
    end
end

%&误差函数
f = 0;
for m = 1:N_sat
    for n = 1:N_sat
        if m == n
            f = f+0;
        else
            f=f+(T_0(n)-(t_inter(m,n)-t_N0(1,m))-t_inter_estimate(m,n))^2;
        end
    end
    f=f+sum((P_0(m,:)-P_0(m,:)).^2);
end

%%各个自变量的倒数
f_t = zeros(1,N_sat);
f_P = zeros(1,3*N_sat);
for n = 1:N_sat
    for m = 1:N_sat
        if m == n
            f_t(1,n) = f_t(1,n)+0;
        else
            f_t(1,n) =  f_t(1,n)+2*(T_0(n)-(t_inter(m,n)-t_N0(1,m))-t_inter_estimate(m,n));
        end
    end
end

for m = 1:N_sat
    for i = 1:3
        for n = 1:N_sat
            if m == n
                f_P(1,3*m-(3-i)) = f_P(1,3*m-(3-i))+0;
            else
                f_P(1,3*m-(3-i)) = f_P(1,3*m-(3-i))+...
                2*(T_0(n)-(t_inter(m,n)-t_N0(1,m))-t_inter_estimate(m,n))*(-1)*((P_0(m,i)-P_0(n,i))/t_inter(m,n)-(P_0(m,i)-P_base(i))/t_N0(1,m))+...
                2*(T_0(m)-(t_inter(n,m)-t_N0(1,n))-t_inter_estimate(n,m))*(-1)*((P_0(m,i)-P_0(n,i))/t_inter(n,m))+2*(P_0(m,i)-P_0(m,i));
            end
        end
    end
end
dj = [f_t,f_P];
end

