 function LADRC1()
%{
�����ܣ�
1��ʹ��m��������LADRC�ĳ���
2��һ�׶���
%}
    %% �����б�
    clear,clc,close all
    global beta_Wo beta_Wc z0 b0 N h
    t1=0;
    t2=5; %����ʱ��
    h=0.001; %��ɢ����
    num=abs(t1-t2)/h ; %����ʱ��ڵ�����
    t=t1: h: t2; %����ʱ��ڵ�    
    N=1; %ϵͳ����
    y0=0; %ϵͳ��ֵ
    z0=zeros(1, N+1);
%     ref=sin(10*t);
    ref=5*ones(1,num+1);
    b0=14;
%     Wo=200;
%     Wc=50;
    Wo=100;
    Wc=20;
   
    
    %% LESOϵ��
    beta_Wo=zeros(N+1,1);
    for i=1 : N+1
        beta_Wo( i )=nchoosek(N+1, i)*( Wo^i );
    end
    
    %% LCϵ��
    beta_Wc=zeros(N,1);
    for i=1 : N
        beta_Wc( i )=nchoosek(N, N+1-i)*( Wc^(N+1-i) );
    end 
    
    %% ���濪ʼ
    y=zeros(num+1, 1); %ϵͳ���ֵԤ��洢�ռ�
    y(1)=y0;
    u=zeros(num+1,1); %ϵͳ������Ԥ��洢�ռ�
    for j=1:num
        
        z=LESO(y(j), u(j), z0);
        z0=z; %�洢��ǰ�ڵ��z����һ�ε���
        u(j+1)=LC(z, ref(j));
        y(j+1)=System(u(j), y(j) );
 
    
    end
    
    %% ������Ӧ����
    figure
    plot(t, y, t ,ref, 'linewidth', 2)
    legend('y', 'ref')
    figure
    plot(t, u, 'linewidth', 2)
    legend('u')
    
end

%% ��������״̬�۲���LESO
function z=LESO(y, u, z0)
    global beta_Wo h  N b0
    % z0��ʾ��һ��ʱ�䲽�����������z
    e=y-z0(1);
    z=zeros(N+1, 1) ;
    for i=1: N-1
        z(i)=z0(i)+h*( beta_Wo(i)*e+z0(i+1) ) ; %LESO΢�ַ�����   
    end
    
    z(N)=z0(N)+h*( beta_Wo(N)*e +z0(N+1)+b0*u );
    z(N+1)=z0(N+1)+h*( beta_Wo(N+1)*e );
    
end

%% ���Է���������LC
function u=LC(z, ref)
    global beta_Wc b0 N
    
    e=z(1)- ref;
    u0=e*beta_Wc(1)+z(end);
    for i=2 : N
        u0=u0+beta_Wc( i )*z(i) ;
    end
   
    u= -u0/b0 ;
    
end

% һ��ϵͳ΢�ַ�������

function y=System(u, yk)
    global h
    y=yk+h*(-10*yk+14*u );

end