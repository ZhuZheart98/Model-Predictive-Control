clear
num = [1];
den = conv([1 1],[1,2]);
sys = tf(num,den);                  %控制对象模型
% figure
% impulse(sys)
%{
===================================
初始化
===================================
%}
T = 0.01;                          %采样时间
N = 1000;                          %采样时域 
P = 100;                          %优化时域
M = 50;                          %控制时域

%内部模型
g = impulse(sys,T:T:N*T);  

%G_1矩阵
G_1 = zeros(P,M);
for i = 1:M-1
    for j = 1:P
        if j-i+1>0
            G_1(j,i) = g(j-i+1);
        end
    end
end
for j = M:P
    G_1(j,M) = G_1(j-1,M-1)+G_1(j-1,M);
end
%{
===================================
第二种构建G_1,利用矩阵
g_1 = g(1:P-M+1);
A = zeros(P-M+1,P-M+1);
for i = 1:P-M+1
    for j = 1:P-M+1
        if i-j+1>0
            A(i,j)=1;
        end
    end
end
A = A*g_1;
for i = M:P
    G_1(i,M) = A(i-M+1);
end
===================================
%}

%G_2矩阵
G_2 = zeros(P,N-1);
for i = 1:P
    for j = 1:N-1
        if i+j<=N
            G_2(i,j) = g(i+j);
        end
    end
end

%误差全矩阵Q，控制权矩阵R
Q = eye(P); R = eye(M); 
C = zeros(1,M);C(1) = 1;         %取首元素的运算
d_1 = inv(G_1'*Q*G_1+R)*G_1'*Q;  %控制向量
d_2 = C*d_1;
h = ones(P,1);                   %校正向量h

%移位阵S
S = zeros(N-1,N-1);
for i = 2:N-1
    S(i,i-1) = 1;
end

%Y_r参考轨迹
Y_r = zeros(P,1);
tau = 0.01;
alpha = exp(-T/tau);

Y_M = zeros(P,1);                %预测值
y = zeros(N,1);                  %实际值
u_1 = zeros(M,1);                %u1，k时刻及k时刻以后的的控制律
u_2 = zeros(N-1,1);              %u2，k时刻以前的控制律
x10 = 0; x20 = 0;    
aa = exp(-T);
bb = exp(-T * 2); 


%{
===================================
实时控制
===================================
%}
for k = 1:N
    %Y_r参考轨迹
    for i = 1:P
        Y_r(i) = alpha^i*y(k)+(1-alpha^i)*5;  %c=y(k)
    end
    
    e(k) = y(k)-Y_M(1);                %误差
    Y_P = Y_M+h.*e(k);                 %预测值校正
    u_1 = d_1*(Y_r-G_2*u_2-h.*e(k));   %最优控制律
    u(k) = C*u_1;
    
    %预测值
    Y_M = G_1*u_1+G_2*u_2;

    %实际仿真模型
    x11 = aa*x10+(1-aa)*u(k)/T;  
    x21 = bb*x20+(1-bb)*x11/2;  
    x10 = x11; x20 = x21;
    y(k+1) = x21;
    
    u_2 = S*u_2;
    u_2(1) = u(k);
end

t = T*(1:N);
subplot(211);
%figure
plot(T*(1:N+1),y);
title('MAC')
grid on
xlabel('Time(s)');
ylabel('Output')
legend('y');

subplot(212);
%figure
plot(t,u);
grid on
xlabel('Time(s)');
legend('u');