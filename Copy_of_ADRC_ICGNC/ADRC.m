clc
clear
%% ADRC参数初始化
h=1/1e3; alpha=1; xd=2; h1=alpha*h; % TD
alpha1=0.5; alpha2=0.25; del=0; b=1; beta01=0; beta02=0; beta03=0; % ESO
beta1=0; beta2=0; gam=1; c=1; h2=gam*h; % NLSEF
%% 跟踪微分器TD
e=xd1(kkk)-xd;
xd1(kkk+1)=xd1(kkk)+h*xd2(kkk);
xd2(kkk+1)=xd2(kkk)+h*fhan(e,xd2(kkk),r,h1);
%% 扩张状态观测器ESO
e=z1(kkk)-y(kkk);
z1(kkk+1)=z1(kkk)+h*(z2(kkk)-beta01*e);
z2(kkk+1)=z2(kkk)+h*(z3(kkk)-beta02*fal(e,alpha1,del)+b*u(kkk));
z3(kkk+1)=z3(kkk)+h*(-beta03*fal(e,alpha2,del));
%% 非线性状态误差反馈律NLSEF
e1=xd1(kkk)-z1(kkk);
e2=xd2(kkk)-z2(kkk);
% u0(kkk)=beta1*fal(e1,alpha1,del)+beta2*fal(e2,alpha2,del);
u0(kkk)=-fhan(e1,c*e2,r,h2);
u(kkk)=(u0(kkk)-z3(kkk))/b;

function fh=fhan(x1, x2, r, h) % fhan函数
    d=r*h; d0=d*h; y=x1+h*x2; a0=sqrt(d^2+8*r*abs(y));
    if(abs(y)>d0), a=x2+sign(y)*(a0-d)/2;
    else, a=x2+y/h;
    end
    if(abs(a)>d), fh=-r*a*sign(a);
    else, fh=-r*a/d;
    end
end
function fe=fal(e, alpha, del) % fal函数描述
    if(abs(e)<=del), fe=e/(del^(1-del));
    else, fe=sign(e)*abs(e)^alpha;
    end
end
