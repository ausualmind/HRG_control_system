function LADRC2()
%{
程序功能：
1、使用m语言描述LADRC的程序
2、二阶对象
%}
    %% 参数列表
    clear,clc,close all
    global beta_Wo beta_Wc z0 b0 N h
    t1=0;
    t2=5; %仿真时间
    h=0.001; %离散步长
    num=abs(t1-t2)/h ; %仿真时间节点数量
    t=t1: h: t2; %仿真时间节点
    
    N=2; %系统阶数
    y0=0; %系统初值
%     dy0=2; %导数初值
    z0=zeros(1, N+1);
%     ref=sin(10*t);
    ref=ones(1,num+1);
    b0=0.04; %和系统b近似
%     Wo=200;
%     Wc=50;
    Wo=100;
    Wc=20;
   
    
    %% LESO系数
    beta_Wo=zeros(N+1,1);
    for i=1 : N+1
        beta_Wo( i )=nchoosek(N+1, i)*( Wo^i );
    end
    
    %% LC系数
    beta_Wc=zeros(N,1);
    for i=1 : N
        beta_Wc( i )=nchoosek(N, N+1-i)*( Wc^(N+1-i) );
    end 
    
    %% 仿真开始
    y=zeros(num+1, 1); %系统输出值预设存储空间
%     dy=zeros(num+1, 1); %系统输出的导数预设存储空间
    y(1)=y0;
    y(2)=0;
    u=zeros(num+1,1); %系统控制量预设存储空间
    for j=1:num-1
        
        z=LESO(y(j), u(j), z0);
        z0=z; %存储当前节点的z，下一次迭代
        u(j+1)=LC(z, ref(j));
        y(j+2)=System(u(j+1), y(j), y(j+1));
    
    end
    
    %% 绘制响应曲线
    figure
    plot(t, y,'--',  t ,ref, 'linewidth', 2)
    legend('y', 'ref')
    figure
    plot(t, u, 'linewidth', 2)
    legend('u')
    figure 
    plot(t, y)

end

%% 线性扩张状态观测器LESO
function z=LESO(y, u, z0)
    global beta_Wo h  N b0
    % z0表示上一个时间步长计算出来的z
    e=y-z0(1);
    z=zeros(N+1, 1) ;
    for i=1: N-1
        z(i)=z0(i)+h*( beta_Wo(i)*e+z0(i+1) ) ; %LESO微分方描述   
    end
    
    z(N)=z0(N)+h*( beta_Wo(N)*e +z0(N+1)+b0*u );
    z(N+1)=z0(N+1)+h*( beta_Wo(N+1)*e );
    
end

%% 线性反馈控制器LC
function u=LC(z, ref)
    global beta_Wc b0 N
    
    e=z(1)- ref;
    u0=e*beta_Wc(1)+z(end);
    for i=2 : N
        u0=u0+beta_Wc( i )*z(i) ;
    end
   
    u= -u0/b0 ;
    
end

%% 一阶系统微分方程描述

% function y=System(u, yk)
%     global h
%     y=yk+h*(-10*yk+14*u );
% 
% end
%% 二阶系统微分方程描述

function y=System(u, yk, yk1)
    global h
 
    y=2*yk1-yk+1/0.9*h*(-yk1+yk)+0.033/0.9*u*h^2;

end