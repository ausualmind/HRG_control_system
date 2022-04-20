function C_PID_01()
%{
�����ܣ�
1��ʹ��m��������PID����ϵͳ
2��

%}
clear,clc,close all
    global h 
%% PID��������������
    kp=9;
    ki=0.2;
    kd=1;
    
%% �����������
    t1=0;
    t2=10;
    h=0.01; %ʱ�䲽��
    num=abs(t1-t2)/h;
    t=t1:h:t2; %����ʱ��ڵ�
    y0=0; %ϵͳ��ֵ
    ref=1*ones(num+1, 1); %�ο�����
    y=zeros(num+1, 1); %ϵͳ���Ԥ��ռ�
    y(1)=y0;
    u=zeros(num,1); %�������洢�ռ�Ԥ��
    
%% ����ϵͳ��������
    
%     forward_error=ref(1)-y(1); %��һ�������
    e=zeros(num+1, 1);
    e(1)=ref(1)-y(1);
    sum_error=0;
    
    for j=1:num
        
        e(j+1)=ref(j)-y(j);
        u(j)=PIDController(e(j+1), e(j), sum_error, kp, ki, kd);
        y(j+1)=System( u(j), y(j) )  ;%�������
        
        sum_error=sum_error+ e(j+1);%���Ļ���
         
        
    end
    
%% ����ͼ��
    plot(t, ref, 'b--', 'linewidth', 2)
    hold on
    plot(t, y, 'k-', 'linewidth', 2)
    title('C and PID Controller')
    legend('ref', 'y')
    axis([t1,t2, 0, 1.2])
    
end
    
%%PID�������Ӻ���
function u=PIDController(e ,fe, sum_error , kp, ki, kd)
    global h
    
%     e=ref-y;  %�������
%     de= (e-fe)/h  ;
    de= (e-fe);
    u=kp*e+ki*sum_error+kd*de;
    
    
end

%% ����ϵͳ����ɢ�����ʽ
function yk1=System(u, yk)
%{
���ݺ�����
G=Y/U=1/(s+1)

%}
    global h
     yk1=h*(u-yk)+yk;

end