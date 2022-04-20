function C_PID_01()
%{
程序功能：
1、使用m语言描述PID控制系统
2、

%}
clear,clc,close all
    global h 
%% PID控制器参数设置
    kp=9;
    ki=0.2;
    kd=1;
    
%% 仿真参数设置
    t1=0;
    t2=10;
    h=0.01; %时间步长
    num=abs(t1-t2)/h;
    t=t1:h:t2; %仿真时间节点
    y0=0; %系统初值
    ref=1*ones(num+1, 1); %参考输入
    y=zeros(num+1, 1); %系统输出预设空间
    y(1)=y0;
    u=zeros(num,1); %控制量存储空间预设
    
%% 控制系统迭代过程
    
%     forward_error=ref(1)-y(1); %上一步的误差
    e=zeros(num+1, 1);
    e(1)=ref(1)-y(1);
    sum_error=0;
    
    for j=1:num
        
        e(j+1)=ref(j)-y(j);
        u(j)=PIDController(e(j+1), e(j), sum_error, kp, ki, kd);
        y(j+1)=System( u(j), y(j) )  ;%被控输出
        
        sum_error=sum_error+ e(j+1);%误差的积分
         
        
    end
    
%% 绘制图像
    plot(t, ref, 'b--', 'linewidth', 2)
    hold on
    plot(t, y, 'k-', 'linewidth', 2)
    title('C and PID Controller')
    legend('ref', 'y')
    axis([t1,t2, 0, 1.2])
    
end
    
%%PID控制器子函数
function u=PIDController(e ,fe, sum_error , kp, ki, kd)
    global h
    
%     e=ref-y;  %跟踪误差
%     de= (e-fe)/h  ;
    de= (e-fe);
    u=kp*e+ki*sum_error+kd*de;
    
    
end

%% 被控系统的离散差分形式
function yk1=System(u, yk)
%{
传递函数：
G=Y/U=1/(s+1)

%}
    global h
     yk1=h*(u-yk)+yk;

end