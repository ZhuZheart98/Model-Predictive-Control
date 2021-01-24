clear 
clc
% num = [1];
% den = conv([1 1],[1,2]);
% sys = tf(num,den);                  %控制对象模型

%{
===================================
初始化
===================================
%}
a = cell(1,2);b = cell(1,2);c = cell(1,1);d = cell(1,1);  %控制对象参数
syms k;
k = length(k);

if k>=0 && k<=150
    a = [1 0.9234]; b = [7.2402 0.9485]; c = 1; d = 1;
elseif k>150 && k<=300
    a = [1 0.8981]; b = [9.9901 0.14142]; c = 1; d = 1;
elseif k>300 && k<=450
    a = [1 0.8838]; b = [9.6041 0.34067]; c = 1; d = 1;
elseif k>450 && k<=600
    a = [1 0.9234]; b = [7.2402 0.9485]; c = 1; d = 1;
end


n_a = length(a)-1;
b = [zeros(1,d-1) b];n_b = length(b)-1;  %n_a,n_b为多项式A,B的阶次（因为d！=1，对b添0）
aa = conv(a,[1 -1]);naa = n_a+1;  %aa的阶次
N1 = d;N = 15;NU = 5; %最小、最大优化时域，控制时域
gamma = 1*eye(NU);alpha = 0.11;%控制加权矩阵、输出柔化系数

L = 600;%控制步数
uk = zeros(d+n_b,1); %输入初值：uk(i)表示u（k-i）
duk = zeros(d+n_b,1);%控制增量
yk = zeros(naa,1);%输入初值
w = 10*[ones(L/4,1);-ones(L/4,1);ones(L/4,1);-ones(L/4+d,1)]; %设定值
xi = sqrt(0.01)*randn(L,1); %白噪声序列

[E,F,G] = multidiophantine(aa,b,c,N);
G = G(N1:N,:);
F1 = zeros(N - N1 + 1, NU);F2 = zeros(N - N1 + 1,n_b);
for i = 1 : N - N1 + 1
    for j = 1 : min(i,NU) 
        F1(i,j) = F(i+N1-1, i+N1-1-j+1);
    end
    for j = 1:n_b
        F2(i,j) = F(i+N1-1,i+N1-1+j);
    end
end



%{
===================================
仿真
===================================
%}

for k = 1:L
    if k>=0 && k<=150
        times(k) = k;
        a = [1 0.9234]; b = [7.2402 0.9485]; c = 1; d = 1;
        y(k) = -aa(2:naa+1)*yk+b*duk(1:n_b+1)+xi(k);
        Yk = [y(k);yk(1:n_a)];
        dUk = duk(1:n_b);
        
    elseif k>150 && k<=300
        times(k) = k;
        a = [1 0.8981]; b = [9.9901 0.14142]; c = 1; d = 1;
        y(k) = -aa(2:naa+1)*yk+b*duk(1:n_b+1)+xi(k);
        Yk = [y(k);yk(1:n_a)];
        dUk = duk(1:n_b);
    
    elseif k>300 && k<=450
        times(k) = k;
        a = [1 0.8838]; b = [9.6041 0.34067]; c = 1; d = 1;
        y(k) = -aa(2:naa+1)*yk+b*duk(1:n_b+1)+xi(k);
        Yk = [y(k);yk(1:n_a)];
        dUk = duk(1:n_b);
        
    elseif k>450 && k<=L
        times(k) = k;
        a = [1 0.9234]; b = [7.2402 0.9485]; c = 1; d = 1;
        y(k) = -aa(2:naa+1)*yk+b*duk(1:n_b+1)+xi(k);
        Yk = [y(k);yk(1:n_a)];
        dUk = duk(1:n_b);
    end
    
    %参考轨迹
    yr(k) = y(k);
    for i = 1:N
        yr(k+i) = alpha*yr(k+i-1)+(1-alpha)*w(k+d);
    end
    Yr = [yr(k+N1:k+N)]';
    
    dU = inv(F1'*F1+gamma)*F1'*(Yr-F2*dUk-G*Yk);
    du(k) = dU(1);u(k) = uk(1)+du(k);
    for i = 1+n_b:-1:2
        uk(i) = uk(i-1);
        duk(i) = duk(i-1);
    end
    uk(1) = u(k);
    duk(1) = du(k);
    for i = naa:-1:2
        yk(i) = yk(i-1);
    end
end

subplot(211);
plot(times,w(1:L),'m:',times,y);
xlabel('k');ylabel('w(k)、y(k)');
legend('w(k)','y(k)');
subplot(212);
plot(times,u);
xlabel('k');ylabel('u(k)');

    
















