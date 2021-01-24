clc
close all;
clear all;
clear
num = [1];
den = conv([1 1],[1,2]);
sys = tf(num,den);                  %���ƶ���ģ��
T = 0.01;                          %����ʱ��
sys_2 = c2d(sys,T);
[num1,den1] = tfdata(sys_2,'v');
%w = 5;
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
Q = 1*eye(P); R = 1*eye(M); 
C = zeros(1,M);C(1) = 1;       %ȡ��Ԫ�ص�����
C_2 = zeros(1,N);C_2(1) = 1;       %ȡ��Ԫ�ص�����
d = C*inv(A'*Q*A+R)*A'*Q;      %��������d
h = 1*ones(N,1);                 %У������h
%��λ��S
S = zeros(N,N);
for i = 1:N-1
    S(i,i+1) = 1;
end
S(N,N) = 1;
ONE = ones(P,1);
I_1 = eye(P,P);
I_2 = zeros(P,N-P);
I = [I_1,I_2];
sim('DMC_continus_system')

subplot(2,1,1);
plot(y,'LineWidth',2);
hold on;
plot(w,':r','LineWidth',2);
xlabel('\fontsize{15}k');
ylabel('\fontsize{15}y,w');
legend('���ֵ','�趨ֵ');
grid on;
subplot(2,1,2);
plot(u,'LineWidth',2)
xlabel('\fontsize{15}k');
ylabel('\fontsize{15}u');
grid on;