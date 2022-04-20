clc
clear
tic
hh=waitbar(0,'please wait');
%% 基本参数
m0=1.52961; m1=0.55296; b0=-0.42369; e0=1.05922; e1=0.31777;
gam=0.17; ome0=4964.703945*2*pi; OME=0/180*pi; dxi0=0.0095429;
%% 误差初始化
eps_rou4=6.13e-4/2200; the_rou4=30; % reference
eps_h4=1.18e-10/(0.85e-3); the_h4=40;
eps_R4=1.05e-9/0.015; the_R4=50;
eps_E4=2.16e4/(7.76e10); the_E4=60;
eps_gam4=5.55e-8/0.17; the_gam4=70;
eps_xi4=0.0347; the_xi4=80;
%% 参数整合
arou4=eps_rou4*cosd(the_rou4); brou4=eps_rou4*sind(the_rou4);
ah4=eps_h4*cosd(the_h4); bh4=eps_h4*sind(the_h4);
aR4=eps_R4*cosd(the_R4); bR4=eps_R4*sind(the_R4);
aE4=eps_E4*cosd(the_E4); bE4=eps_E4*sind(the_E4);
agam4=eps_gam4*cosd(the_gam4); bgam4=eps_gam4*sind(the_gam4);
axi4=eps_xi4*cosd(the_xi4); bxi4=eps_xi4*sind(the_xi4);

ke1=3*m1*ah4-4*m1*aR4+m1*aE4-m1*agam4/(1+gam);
ke2=3*m1*bh4-4*m1*bR4+m1*bE4-m1*bgam4/(1+gam);
ke3=arou4+ah4;
ke4=brou4+bh4;
Cpp=(m0*ome0^2*ke1-(m1*(m0*ome0^2-e0*OME^2)+m0*e1*OME^2)*ke3)/...
    (m0^2-m1^2*(ke3^2+ke4^2));
Cpq=(-m1*ome0^2*(ke2*ke3-ke1*ke4))/...
    (m0^2-m1^2*(ke3^2+ke4^2));
Epp=(m1*e1*OME^2*(ke3^2+ke4^2)-m1*ome0^2*(ke1*ke3+ke2*ke4))/...
    (m0^2-m1^2*(ke3^2+ke4^2))+...
    (m1^2*(m0*ome0^2-e0*OME^2)*(ke3^2+ke4^2))/...
    (m0^3);
Epq=(m0*ome0^2*ke2-(m0*e1*OME^2+m1*(m0*ome0^2-e0*OME^2))*ke4)/...
    (m0^2-m1^2*(ke3^2+ke4^2));
Dpp=(dxi0*(m1*axi4+ke1)-dxi0*m1*ke3+4*b0*m1*OME*ke4)/...
    (m0^2-m1^2*(ke3^2+ke4^2));
Dpq=(-dxi0*m1^2/m0*(bxi4*ke3-axi4*ke4)-dxi0*m1/m0*(ke2*ke3-ke1*ke4))/...
    (m0^2-m1^2*(ke3^2+ke4^2))-...
    (4*b0*OME*m1^2*(ke3^2+ke4^2))/...
    (m0^3);
Gpp=(-dxi0*m1^2/m0*(axi4*ke3+bxi4*ke4)-dxi0*m1/m0*(ke1*ke3+ke2*ke4))/...
    (m0^2-m1^2*(ke3^2+ke4^2))+(dxi0*m0)/...
    (m0^2-m1^2*(ke3^2+ke4^2));
Gpq=(-4*b0*OME*m1*ke3+dxi0*(m1*bxi4+ke2-m1*ke4))/...
    (m0^2-m1^2*(ke3^2+ke4^2));
%% 动力学状态方程矩阵
a=-(Dpp+Gpp);
b=-4*b0*OME/m0-(Dpq+Gpq);
c=-(ome0^2-e0*OME^2/m0)-(Cpp+Epp);
d=-(Cpq+Epq);
e=4*b0*OME/m0+(Dpq-Gpq);
f=Dpp-Gpp;
g=Cpq-Epq;
h=-(ome0^2-e0*OME/m0)+(Cpp-Epp);
A=[a b c d
   e f g h
   1 0 0 0
   0 1 0 0];
xx=9.70537e-8; del_phiv=0.5/180*pi;
V1a4=0.01; V1b4=0.01; V2a4=0.01; V2b4=0.01;
r=xx*(1+V1a4); s=xx*(V2b4-2*del_phiv);
u=xx*V1b4; v=xx*(1-V2a4);
B=[r s
   u v
   0 0
   0 0];
C=[0 0 1 0
   0 0 0 1];
%% 求解G(k)
syms s t Ls;   % 求状态转移矩阵 利用拉氏变换，syms为符号函数用来定义数学函数
    I = eye(size(A));
    Ls = inv(s*I - A); % collect 函数为合并同类项
    STM = ilaplace(Ls,s,t); %状态转移矩阵,ilaplace为拉氏反变换函数
    G = double(subs(STM,t,1e-6)); % 符号函数求解
%% 求解H(k)
syms T
    HLs = int(STM,t,0,T);
    H = HLs*B;
    H = real(double(subs(H,T,1e-6))); % 符号函数求解
%% 动力学
TT=20e6; xxx=zeros(4,TT+1); FF=zeros(2,TT+1); phi=zeros(1,TT+1);
x=zeros(1,TT+1); y=zeros(1,TT+1); p=zeros(1,TT+1); q=zeros(1,TT+1);
x(1)=10e-6*ome0; F1=zeros(1,TT+1); F2=zeros(1,TT+1);
pVc=zeros(1,TT/8+1); qVc=zeros(1,TT/8+1);
pVs=zeros(1,TT/8+1); qVs=zeros(1,TT/8+1);
Cp=zeros(1,TT/8+1); Cq=zeros(1,TT/8+1);
Sp=zeros(1,TT/8+1); Sq=zeros(1,TT/8+1);
S_delphi=zeros(1,TT/1e3+1); C_delphi=zeros(1,TT/1e3+1);
fp=zeros(1,TT/1e3+1); fome=zeros(1,TT/1e3+1);
del_ome=zeros(1,TT/1e3+1); 
for iteration=1:1
    for k=1:TT % 滤波器初始化
        kk=idivide(k,int32(8),'ceil');
        kkk=idivide(k,int32(1e3),'ceil');
        if mod((k),1e5)==0
            str=['运行中...',num2str(k/1e6),'秒'];
            waitbar(k/TT,hh,str)
        end
        xxx(:,k)=[x(k) y(k) p(k) q(k)]'; FF(:,k)=[F1(k) F2(k)]';
        x(k+1)=G(1,:)*xxx(:,k) + H(1,:)*FF(:,k);
        y(k+1)=G(2,:)*xxx(:,k) + H(2,:)*FF(:,k);
        p(k+1)=G(3,:)*xxx(:,k) + H(3,:)*FF(:,k);
        q(k+1)=G(4,:)*xxx(:,k) + H(4,:)*FF(:,k);
        %%% 解调滤波
        if mod(k,8)==0
            pVc(kk)=2e5*p(k)*2*cos(phi(k)); qVc(kk)=2e5*q(k)*2*cos(phi(k));%正交
            pVs(kk)=2e5*p(k)*2*sin(phi(k)); qVs(kk)=2e5*q(k)*2*sin(phi(k));%同相
            b1=-1.997238; b2=0.997245;
            a0=0.000101412; a1=-0.00019662; a2=0.000101411;
        end
        if (idivide(k,int32(8))>2)&&(mod(k,8)==0)
            Cp(kk)=-b1*Cp(kk-1)-b2*Cp(kk-2)+a0*pVc(kk)+a1*pVc(kk-1)+a2*pVc(kk-2);
            Cq(kk)=-b1*Cq(kk-1)-b2*Cq(kk-2)+a0*qVc(kk)+a1*qVc(kk-1)+a2*qVc(kk-2);
            Sp(kk)=-b1*Sp(kk-1)-b2*Sp(kk-2)+a0*pVs(kk)+a1*pVs(kk-1)+a2*pVs(kk-2);
            Sq(kk)=-b1*Sq(kk-1)-b2*Sq(kk-2)+a0*qVs(kk)+a1*qVs(kk-1)+a2*qVs(kk-2);
        end
        %%% 多PI/频率跟踪算法
        if (idivide(k,int32(1e3))>1)&&(mod(k,8)==0)
            Ts3=1/1e3;
            S_delphi(kkk)=2*(Cp(kk)*Sp(kk)+Cq(kk)*Sq(kk));
            C_delphi(kkk)=Cp(kk)^2+Cq(kk)^2-Sp(kk)^2-Sq(kk)^2;
            del_ome(kkk)=1/2*jdjs(-S_delphi(kkk),-C_delphi(kkk));
            P1=0.2; I1=1; P2=40; I2=2000;
%             fome(kkk)=fome(kkk-1)+P2*(-del_ome(kkk))+(I2*Ts3-P2)*(-del_ome(kkk-1));
            fp(kkk)=fp(kkk-1)+P1*(2-Sp(kk)*1.122)+(I1*Ts3-P1)*(2-Sp(kk-1)*1.122);
        end
        %%% 调制
        phi(k+1)=phi(k)+(fome(kkk)+ome0)/1e6;
        F1(k+1)=1/xx*fp(kkk)*cos(phi(k+1));
    end
end
delete(hh);
toc
function y = jdjs(S,C)
    y = 0;
    if S >= 0
        if C > 0
            y = atan(S/C);
        end
        if C < 0
            y = (pi + atan(S/C) );
        end
        if C == 0
            y = pi/2;
        end
    else
        if C > 0
            y = ( atan(S/C) );
        end
        if C < 0
            y = ( atan(S/C) - pi);
        end
        if C == 0
            y = -pi/2;
        end
    end
end
