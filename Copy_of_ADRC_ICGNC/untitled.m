clc
clear
tic
%% HRG��ʼ��
GHox=HRG_initia();
TTTT=1/40e3; % HRGʰȡ����ʱ��
%% ���G(k)
syms s t Ls;   % ��״̬ת�ƾ��� �������ϱ任��symsΪ���ź�������������ѧ����
    I = eye(size(GHox.A));
    Ls = inv(s*I - GHox.A); % collect ����Ϊ�ϲ�ͬ����
    STM = ilaplace(Ls,s,t); %״̬ת�ƾ���,ilaplaceΪ���Ϸ��任����
    GHox.GG = double(subs(STM,t,TTTT)); % ���ź������
%% ���H(k)
syms T
    HLs = int(STM,t,0,T);
    GHox.HH = HLs*GHox.B;
    GHox.HH = real(double(subs(GHox.HH,T,TTTT))); % ���ź������
%% PI������ʼ��
P1=40; I1=2000; P2=0.2; I2=1; P3=35; I3=2500; P4=35; I4=2500; % PI����������
PIs=[P1,I1;P2,I2;P3,I3;P4,I4];
%% ADRC������ʼ��
ADRC.h=1/1e3; ADRC.r=100; ADRC.alpha=1; ADRC.h1=ADRC.alpha*ADRC.h; % TD
ADRC.alpha1=0.5; ADRC.alpha2=0.25; ADRC.del=0.01; ADRC.b=1; % ESO
ADRC.gam=1; ADRC.c=1; ADRC.h2=ADRC.gam*ADRC.h; % ��Щ�������ã�
ADRC.flag12=1; % flag12=1��ʾʹ��beta12 / flag12=2��ʾʹ��gam,c,h2

ADRC.alpha3=0.1; ADRC.alpha4=1;
ADRC.beta01=1000; ADRC.beta02=1.094743921412483; ADRC.beta03=0.18585878393097;
ADRC.beta1=0; ADRC.beta2=0; % ����PD���� % NLSEF
% ADRC.beta01=1000; ADRC.beta02=1.094743921412483; ADRC.beta03=0.18585878393097;
% ADRC.beta1=0; ADRC.beta2=0; % ����PD���� % NLSEF
%% PSO������ʼ��
inref.noP=36; inref.nn=1; % PSO paramters
inref.Max_iteration = 200; inref.Vmax=[0.5,0.05];
inref.wmax=0.9; inref.wmin=0.2;
inref.c1=10; inref.c2=10;
% ������ԭ�����ʱ��ע���Ŀǰ����Ѱ��״̬��2�����ĸ�����PSO_w����33/102��
inref.X1min=0; inref.X1max=2; % X1/X2�������½� // PSO��
inref.X2min=0; inref.X2max=0.2; % PSO��
%% PSO/������/HRG������· ���Ϸ���
GHox.err1234=2; % [err_ome,err_ps,err_qc,err_qs] / ѡ����Ƶ��������
preset.ps.xd=2; preset.qc.xd=0; preset.qs.xd=0; % Ԥ��
fobj = @HRG_con;    
con_style=2; % 1��ʾѡ��PI / 2��ʾѡ��ADRC
[Swarm,cg_curve] = PSO_w(inref,fobj,PIs,ADRC,preset,GHox,con_style);
% J = HRG_con(PIs,ADRC,preset,GHox,con_style);

toc
