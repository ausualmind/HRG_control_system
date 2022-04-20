function GHox = HRG_initia() % HRG初始化
%% 基本参数
m0=1.52961; m1=0.55296; b0=-0.42369; e0=1.05922; e1=0.31777;
gam=0.17; GHox.ome0=4964.703945*2*pi; OME=20/180*pi;
dxi0=0.00954297824972570791914581868957; % 阻尼项的值
%% 误差初始化
eps_rou4=6.13e-4/2200; the_rou4=30; % 密度不均匀
eps_h4=1.18e-10/(0.85e-3); the_h4=40; % 厚度不均匀
eps_R4=1.05e-9/0.015; the_R4=50; % 半径不均匀
eps_E4=2.16e4/(7.76e10); the_E4=60; % 弹性模量不均匀
eps_gam4=5.55e-8/0.17; the_gam4=70; % 泊松比不均匀
eps_xi4=0.0347; the_xi4=80; % eps_xi4=0.0347 阻尼不均匀
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
Cpp=(m0*GHox.ome0^2*ke1-(m1*(m0*GHox.ome0^2-e0*OME^2)+m0*e1*OME^2)*ke3)/...
    (m0^2-m1^2*(ke3^2+ke4^2));
Cpq=(-m1*GHox.ome0^2*(ke2*ke3-ke1*ke4))/...
    (m0^2-m1^2*(ke3^2+ke4^2));
Epp=(m1*e1*OME^2*(ke3^2+ke4^2)-m1*GHox.ome0^2*(ke1*ke3+ke2*ke4))/...
    (m0^2-m1^2*(ke3^2+ke4^2))+...
    (m1^2*(m0*GHox.ome0^2-e0*OME^2)*(ke3^2+ke4^2))/...
    (m0^3);
Epq=(m0*GHox.ome0^2*ke2-(m0*e1*OME^2+m1*(m0*GHox.ome0^2-e0*OME^2))*ke4)/...
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
c=-(GHox.ome0^2-e0*OME^2/m0)-(Cpp+Epp);
d=-(Cpq+Epq);
e=4*b0*OME/m0+(Dpq-Gpq);
f=Dpp-Gpp;
g=Cpq-Epq;
h=-(GHox.ome0^2-e0*OME/m0)+(Cpp-Epp);
GHox.A=[a b c d
        e f g h
        1 0 0 0
        0 1 0 0];
GHox.xx=9.70537e-8/(m0^2-(m1*ke3)^2-(m1*ke4)^2);
del_phiv=0/180*pi; V1a4=0; V1b4=0; V2a4=0; V2b4=0;
r=GHox.xx*( (m0-m1*ke3)*(1+V1a4)-m1*ke4*V1b4 );
s=GHox.xx*( (m0-m1*ke3)*(V2b4-2*del_phiv)-m1*ke4*(1-V2a4) );
u=GHox.xx*( (m0+m1*ke3)*V1b4-m1*ke4*(1+V1a4) );
v=GHox.xx*( (m0+m1*ke3)*(1-V2a4)-m1*ke4*(V2b4-2*del_phiv) );
GHox.B=[r s
        u v
        0 0
        0 0];
GHox.C=[0 0 1 0
        0 0 0 1];
end
