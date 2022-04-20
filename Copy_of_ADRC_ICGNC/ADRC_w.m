function [u,xd12,z123] = ADRC_w(ADRC,yy,preset,xd12,z123,u) % HRG控制电路
%% ADRC参数初始化
h=ADRC.h; r=ADRC.r; alpha=ADRC.alpha; xd=preset.xd; h1=ADRC.h1; % TD
alpha1=ADRC.alpha1; alpha2=ADRC.alpha2; del=ADRC.del; b=ADRC.b; % ESO
beta01=ADRC.beta01; beta02=ADRC.beta02; beta03=ADRC.beta03;
beta1=ADRC.beta1; beta2=ADRC.beta2; % 类似PD参数 % NLSEF
gam=ADRC.gam; c=ADRC.c; h2=gam*h;
%% 跟踪微分器TD
xd1=xd12(1,:); xd2=xd12(2,:);
e=-(xd-xd1(1));
xd1(2)=xd1(1)+h*xd2(1);
% xd2(2)=xd2(1)+h*fhan(e,xd2(1),r,h1); % 非线性
xd2(2)=xd2(1)+h*( -r^2*e-2*r*xd2(1) ); % 线性
%% 扩张状态观测器ESO
z1=z123(1,:); z2=z123(2,:); z3=z123(3,:);
e=z1(1)-yy;
% alpha1=1; alpha2=1; % 线性
z1(2)=z1(1)+h*(z2(1)-beta01*e);
z2(2)=z2(1)+h*(z3(1)-beta02*fal(e,alpha1,del)+b*u(1));
z3(2)=z3(1)+h*(-beta03*fal(e,alpha2,del));
%% 非线性状态误差反馈律NLSEF 算法选择
e1=xd1(2)-z1(2); e2=xd2(2)-z2(2);
if ADRC.flag12==1, u0=beta1*fal(e1,ADRC.alpha3,del)+beta2*fal(e2,ADRC.alpha4,del);
else
    if ADRC.flag12==2, u0=-fhan(e1,c*e2,r,h2); else, u0=0; end
end
u(2)=(u0-z3(2))/b;
%% 更新输出信息
xd12=[xd1;xd2];
z123=[z1;z2;z3];
end
%% fhan函数/fal函数
function fh=fhan(x1, x2, r, h) % fhan函数
    d=r*h; d0=d*h; y=x1+h*x2; a0=sqrt(d^2+8*r*abs(y));
    if(abs(y)>d0), a=x2+(a0-d)/2;
    else, a=x2+y/h;
    end
    if(abs(a)>d), fh=-r*sign(a);
    else, fh=-r*a/d;
    end
end
function fe=fal(e, alpha, del) % fal函数描述
    if(abs(e)<=del), fe=e/(del^(1-alpha));
    else, fe=sign(e)*abs(e)^alpha;
    end
end
