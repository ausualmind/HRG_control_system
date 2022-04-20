clc
clear
tic
%% HRG初始化
GHox=HRG_initia();
TTTT=1/40e3; % HRG拾取采样时间
%% 求解G(k)
syms s t Ls;   % 求状态转移矩阵 利用拉氏变换，syms为符号函数用来定义数学函数
    I = eye(size(GHox.A));
    Ls = inv(s*I - GHox.A); % collect 函数为合并同类项
    STM = ilaplace(Ls,s,t); %状态转移矩阵,ilaplace为拉氏反变换函数
    GHox.GG = double(subs(STM,t,TTTT)); % 符号函数求解
%% 求解H(k)
syms T
    HLs = int(STM,t,0,T);
    GHox.HH = HLs*GHox.B;
    GHox.HH = real(double(subs(GHox.HH,T,TTTT))); % 符号函数求解
%% PI参数初始化
P1=40; I1=2000; P2=0.2; I2=1; P3=35; I3=2500; P4=35; I4=2500; % PI控制器参数
PIs=[P1,I1;P2,I2;P3,I3;P4,I4];
%% ADRC参数初始化
ADRC.h=1/1e3; ADRC.r=100; ADRC.alpha=1; ADRC.h1=ADRC.alpha*ADRC.h; % TD
ADRC.alpha1=0.5; ADRC.alpha2=0.25; ADRC.del=0.01; ADRC.b=1; % ESO
ADRC.gam=1; ADRC.c=1; ADRC.h2=ADRC.gam*ADRC.h; % 这些参数不用！
ADRC.flag12=1; % flag12=1表示使用beta12 / flag12=2表示使用gam,c,h2

ADRC.alpha3=0.1; ADRC.alpha4=1;
ADRC.beta01=1000; ADRC.beta02=1.094743921412483; ADRC.beta03=0.18585878393097;
ADRC.beta1=0; ADRC.beta2=0; % 类似PD参数 % NLSEF
% ADRC.beta01=1000; ADRC.beta02=1.094743921412483; ADRC.beta03=0.18585878393097;
% ADRC.beta1=0; ADRC.beta2=0; % 类似PD参数 % NLSEF
%% PSO参数初始化
inref.noP=36; inref.nn=1; % PSO paramters
inref.Max_iteration = 200; inref.Vmax=[0.5,0.05];
inref.wmax=0.9; inref.wmin=0.2;
inref.c1=10; inref.c2=10;
% 分离性原则求参时，注意对目前处于寻优状态的2参数的更换“PSO_w”中33/102行
inref.X1min=0; inref.X1max=2; % X1/X2参数上下界 // PSO用
inref.X2min=0; inref.X2max=0.2; % PSO用
%% PSO/控制器/HRG开环电路 联合仿真
GHox.err1234=2; % [err_ome,err_ps,err_qc,err_qs] / 选择控制的误差类型
preset.ps.xd=2; preset.qc.xd=0; preset.qs.xd=0; % 预设
fobj = @HRG_con;    
con_style=2; % 1表示选择PI / 2表示选择ADRC
[Swarm,cg_curve] = PSO_w(inref,fobj,PIs,ADRC,preset,GHox,con_style);
% J = HRG_con(PIs,ADRC,preset,GHox,con_style);

toc
