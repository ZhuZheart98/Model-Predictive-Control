clear
num = [1];
den = conv([1 1],[1,2]);
sys = tf(num,den);                  %控制对象模型
figure
step(sys)
T = 0.01;                          %采样时间
w = 10;
N = 700;                          %采样时域  7s稳定
P = 100;                          %优化时域
M = 50;                          %控制时域

%{
===================================
矩阵确定
===================================
%}

%模型参数a
a = step(sys,T:T:N*T);         
%A矩阵
A=zeros(P,M);
for i=1:P
    for j=1:M
        if i-j+1>0
            A(i,j)=a(i-j+1);
        end
    end
end
%误差全矩阵Q，控制权矩阵R
Q = eye(P); R = eye(M); 
C = zeros(1,M);C(1) = 1;       %取首元素的运算
d = C*inv(A'*Q*A+R)*A'*Q;      %控制向量d
h = ones(N,1);                 %校正向量h
%移位阵S
S = zeros(N,N);
for i = 1:N-1
    S(i,i+1) = 1;
end
S(N,N) = 1;
%期待阵W
W = zeros(P,1);
for i = 1:P
    W(i) = w;
end
%单位阵
ONE = ones(P,1);
I_1 = ones(P,P);
I_2 = zeros(P,N-P);
I = [I_1,I_2];
%{
===================================
初始化
===================================
%}
Y = zeros(N,1);     %预测值
y = zeros(N,1);     %实际值
u = zeros(N,1);     %控制量
x10 = 0; x20 = 0;    
aa = exp(-T);
bb = exp(-T * 2); 

%{
===================================
实时控制
===================================
%}

for k = 1:N
    e = y(k)-Y(1);                %误差
    Y = Y+h.*e;                   %预测值校正
    Y = S*Y;                      %移位
    Delta = d*[W-Y(1:P)];         %控制增量delta
    Delta_0(k) = Delta;   
    if k-1 <= 0
        u(k) = Delta;
    else
        u(k) = u(k-1)+Delta;
    end
    %预测值
    Y = Y+a.*Delta;

    %实际仿真模型
    x11 = aa*x10+(1-aa)*u(k);  
    x21 = bb*x20+(1-bb)*x11/2;  
    x10 = x11; x20 = x21;
    y(k+1) = x21;
end
t = T*(1:N);
%subplot(211);
figure
plot(T*(1:N+1),y);
title('DMC')
grid on
xlabel('Time(s)');
ylabel('Output')
legend('y');
%subplot(212);
figure
plot(t,u,t,Delta_0);
grid on
xlabel('Time(s)');
legend('u','\Deltau');