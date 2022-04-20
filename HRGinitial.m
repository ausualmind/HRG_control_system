clc
clear
P3 = 40;    I3 = 2000;  %���໷PI����������
% P4 = 7000000;  I4 = 200000;  
P4 = 35;  I4 = 2500;  %��ƽ��PI����������
P1 = 0.2;    I1 = 1;  %��ֵPI����������
P2 = 35;  I2 = 2500;  %����PI����������
n = 2;  %������
K = 1.122;  %�˲�����������
%% �������
% global m0 m1 b0 e0 e1 omega k E gamma h Q ava_damp_x ava_damp ...
%     delta_damp4 ava_rou delta_rou4 ava_R delta_R4
m0=1.52961; m1=0.55296; b0=-0.42369; e0=1.05922; e1=0.31777;
omega = 4964.703945*2*pi; k = 0.277; %����Ƶ�ʺͽ�������
E=7.76e10; gamma=0.17; h=0.85e-3; %����ģ�������ɱȺͺ��
OMEGA = 0/180*pi; azi_theta=0 ; %���ٶ�����ֵ��פ����λ��
Q=5e6; %Ʒ������ƽ��ֵ
%% ���᲻����
ava_damp_x=0.0015188;
ava_damp=0.0031194; delta_damp4 = 1.08195e-4;  %�������Ϊ3.46846%
theta_damp4 = 0;  %ʱ�䳣��Ϊtau_x�ġ���������ᡱ��X��֮��ļн�
a_damp4=delta_damp4*cos(4*theta_damp4);
b_damp4=delta_damp4*sin(4*theta_damp4);
%% �ܶȲ�����
ava_rou=2200; delta_rou4=6.13e-4; %�ܶ����Ϊ0.27864ppm
theta_rou4=0;
a_rou4=delta_rou4*cos(4*theta_rou4);
b_rou4=delta_rou4*sin(4*theta_rou4);
%% �뾶������
ava_R=0.015; delta_R4=1.05e-9; %�뾶���Ϊ0.07ppm
theta_R4=0;
a_R4=delta_R4*cos(4*theta_R4);
b_R4=delta_R4*sin(4*theta_R4);
%% ״̬���̹���
A11=1/(-(m0+m1*a_rou4)*(m0-m1*a_rou4)+(m1*b_rou4)^2)*...
    ((m0-m1*a_rou4)*ava_damp_x*(1+m1*a_damp4/m0-4*m1*a_R4/m0)-...
    m1*b_rou4*(-4*OMEGA*b0+ava_damp_x*(m1*b_damp4/m0-4*m1*b_R4/m0)));
A12=1/(-(m0+m1*a_rou4)*(m0-m1*a_rou4)+(m1*b_rou4)^2)*...
    ((m0-m1*a_rou4)*(4*OMEGA*b0+ava_damp_x*(m1*b_damp4/m0-4*m1*b_R4/m0))-...
    m1*b_rou4*ava_damp_x*(1-m1*a_damp4/m0+4*m1*a_R4/m0));
A13=1/(-(m0+m1*a_rou4)*(m0-m1*a_rou4)+(m1*b_rou4)^2)*...
    ((m0-m1*a_rou4)*((m0-4*m1*a_R4)*omega^2-(e0+e1*a_rou4)*OMEGA^2)-...
    m1*b_rou4*(-4*m1*b_R4*omega^2-e1*b_rou4*OMEGA^2));
A14=1/(-(m0+m1*a_rou4)*(m0-m1*a_rou4)+(m1*b_rou4)^2)*...
    ((m0-m1*a_rou4)*(-4*m1*b_R4*omega^2-e1*b_rou4*OMEGA^2)-...
    m1*b_rou4*((m0+4*m1*a_R4)*omega^2-(e0-e1*a_rou4)*OMEGA^2));
A21=1/(-(m1*b_rou4)^2+(m0+m1*a_rou4)*(m0-m1*a_rou4))*...
    (m1*b_rou4*ava_damp_x*(1+m1*a_damp4/m0-4*m1*a_R4/m0)-...
    (m0+m1*a_rou4)*(-4*OMEGA*b0+ava_damp_x*(m1*b_damp4/m0-4*m1*b_R4/m0)));
A22=1/(-(m1*b_rou4)^2+(m0+m1*a_rou4)*(m0-m1*a_rou4))*...
    (m1*b_rou4*(4*OMEGA*b0+ava_damp_x*(m1*b_damp4/m0-4*m1*b_R4/m0))-...
    (m0+m1*a_rou4)*ava_damp_x*(1-m1*a_damp4/m0+4*m1*a_R4/m0));
A23=1/(-(m1*b_rou4)^2+(m0+m1*a_rou4)*(m0-m1*a_rou4))*...
    (m1*b_rou4*((m0-4*m1*a_R4)*omega^2-(e0+e1*a_rou4)*OMEGA^2)-...
    (m0+m1*a_rou4)*(-4*m1*b_R4*omega^2-e1*b_rou4*OMEGA^2));
A24=1/(-(m1*b_rou4)^2+(m0+m1*a_rou4)*(m0-m1*a_rou4))*...
    (m1*b_rou4*(-4*m1*b_R4*omega^2-e1*b_rou4*OMEGA^2)-...
    (m0+m1*a_rou4)*((m0+4*m1*a_R4)*omega^2-(e0-e1*a_rou4)*OMEGA^2));
A31=1; A32=0; A33=0; A34=0; A41=0; A42=1; A43=0; A44=0;
B11=1/(-(m0+m1*a_rou4)*(m0-m1*a_rou4)+(m1*b_rou4)^2)*(-1)*(m0-m1*a_rou4)*(m0+m1*a_rou4);
B12=1/(-(m0+m1*a_rou4)*(m0-m1*a_rou4)+(m1*b_rou4)^2)*m1*b_rou4*(m0-m1*a_rou4);
B21=1/(-(m1*b_rou4)^2+(m0+m1*a_rou4)*(m0-m1*a_rou4))*(-1)*m1*b_rou4*(m0+m1*a_rou4);
B22=1/(-(m1*b_rou4)^2+(m0+m1*a_rou4)*(m0-m1*a_rou4))*(m0+m1*a_rou4)*(m0-m1*a_rou4);
B31=0; B32=0; B41=0; B42=0;
%%
t = 0;	%��ʼʱ��
phi0 =  0;	%�ʵ��˶��ĳ�ʼ��λ��
theta = 0;	%������
a=0;% a = 10*10^(-6)*0;	%��Բ�ĳ�����
pq = 1/sqrt(2)*10*10^-6*0;
ka = 1;
q = 0;	%��Բ�Ķ̰���
x0  = [a*cos(omega*t+phi0)*cos(2*theta)-q*sin(omega*t+phi0)*sin(2*theta)
       a*cos(omega*t+phi0)*sin(2*theta)+q*sin(omega*t+phi0)*cos(2*theta)
       -omega*(a*sin(omega*t+phi0)*cos(2*theta)+q*cos(omega*t+phi0)*sin(2*theta))
       -omega*(a*sin(omega*t+phi0)*sin(2*theta)-q*cos(omega*t+phi0)*cos(2*theta))];
pq0  = [pq
       -pq
       0
       0];
%%
A = [A11 A12 A13 A14
     A21 A22 A23 A24
     A31 A32 A33 A34 
     A41 A42 A43 A44];
B = [B11 B12
     B21 B22
     B31 B32
     B41 B42];
C = [0 0 1 0
     0 0 0 1];
Qc = ctrb(A,B);
nc = rank(Qc);
L=length(A);
if nc==L
    str1 = 'ϵͳ״̬��ȫ�ܿ�'
else
    str1 = 'ϵͳ״̬����ȫ�ܿ�'
end
Qo = obsv(A,C);
no = rank(Qo);
if no==L
    str2 = 'ϵͳ״̬��ȫ�ܹ۲�'
else
    str2 = 'ϵͳ״̬����ȫ�ܹ۲�'
end
