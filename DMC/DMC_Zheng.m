clear
num = [1];
den = conv([1 1],[1,2]);
sys = tf(num,den);                  %���ƶ���ģ��
figure
step(sys)
T = 0.01;                          %����ʱ��
w = 10;
N = 700;                          %����ʱ��  7s�ȶ�
P = 100;                          %�Ż�ʱ��
M = 50;                          %����ʱ��

%{
===================================
����ȷ��
===================================
%}

%ģ�Ͳ���a
a = step(sys,T:T:N*T);         
%A����
A=zeros(P,M);
for i=1:P
    for j=1:M
        if i-j+1>0
            A(i,j)=a(i-j+1);
        end
    end
end
%���ȫ����Q������Ȩ����R
Q = eye(P); R = eye(M); 
C = zeros(1,M);C(1) = 1;       %ȡ��Ԫ�ص�����
d = C*inv(A'*Q*A+R)*A'*Q;      %��������d
h = ones(N,1);                 %У������h
%��λ��S
S = zeros(N,N);
for i = 1:N-1
    S(i,i+1) = 1;
end
S(N,N) = 1;
%�ڴ���W
W = zeros(P,1);
for i = 1:P
    W(i) = w;
end
%��λ��
ONE = ones(P,1);
I_1 = ones(P,P);
I_2 = zeros(P,N-P);
I = [I_1,I_2];
%{
===================================
��ʼ��
===================================
%}
Y = zeros(N,1);     %Ԥ��ֵ
y = zeros(N,1);     %ʵ��ֵ
u = zeros(N,1);     %������
x10 = 0; x20 = 0;    
aa = exp(-T);
bb = exp(-T * 2); 

%{
===================================
ʵʱ����
===================================
%}

for k = 1:N
    e = y(k)-Y(1);                %���
    Y = Y+h.*e;                   %Ԥ��ֵУ��
    Y = S*Y;                      %��λ
    Delta = d*[W-Y(1:P)];         %��������delta
    Delta_0(k) = Delta;   
    if k-1 <= 0
        u(k) = Delta;
    else
        u(k) = u(k-1)+Delta;
    end
    %Ԥ��ֵ
    Y = Y+a.*Delta;

    %ʵ�ʷ���ģ��
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