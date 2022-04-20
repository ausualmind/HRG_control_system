clc
clear
tic
h=waitbar(0,'please wait');
% for i=1:1000
%     %computation here%
%     str=['运行中...',num2str(i/1000*100),'%'];
%     waitbar(i/1000,h,str)
% end
sim_T=2;
t=0:1e-6:sim_T;
r_x=zeros(sim_T*1e6+1,1);r_y=zeros(sim_T*1e6+1,1);
m0=1.52961;d_ksi0=0.00954298;ome0=4964.703945*2*pi;
g_xy=2*m0*exp(d_ksi0/m0)*sinh(t*sqrt(d_ksi0^2-4*ome0^2*...
m0^2)/2/m0)/sqrt(d_ksi0^2-4*ome0^2*m0^2);% 传函
xx=zeros(sim_T*1e6+1,1);xx(1)=2;
yy=zeros(sim_T*1e6+1,1); % detect dual init
lpf_amp=zeros(sim_T/8*1e6+1,1);lpf_LOOP=zeros(sim_T/8*1e6+1,1);
lpf_yc=zeros(sim_T/8*1e6+1,1);lpf_ys=zeros(sim_T/8*1e6+1,1); % 解调init
g_LPF1=zeros(sim_T/8*1e6+1,1); % LPF
P1=0.2*2e5;I1=1*2e5;P2=35*2e5;I2=2500*2e5;P3=35*2e5;I3=2500*2e5;A0=2;
f_xc=zeros(sim_T*1e3+1,1);f_yc=zeros(sim_T*1e3+1,1);f_ys=zeros(sim_T*1e3+1,1);
PD=zeros(sim_T*1e3+1,1);delta_ome_r=zeros(sim_T*1e3+1,1); % 多PI
P4=40;I4=2000; % 频率跟踪算法
f_xq=zeros(sim_T*1e6+1,1);f_yq=zeros(sim_T*1e6+1,1); % 调制
for iteration=1:1
    for ii=2:sim_T/8*1e6+1 % 滤波器初始化
        g_LPF1(ii)=-1.396296684*10^(-7)*(0.9986190000+0.002256731929i)^ii-...
            0.001374260186i*(0.9986190000 + 0.002256731929i)^ii -...
            1.396296684*10^(-7)*(0.9986190000 - 0.002256731929i)^ii +...
            0.001374260186i*(0.9986190000 - 0.002256731929i)^ii;
    end
    g_LPF1(1)=0.0001016912593;
    P1=0.2;I1=1;P2=35;I2=2500;P3=35;I3=2500; % 更新电路参数
    for kk=1:sim_T*1e6+1
        if mod((kk),1e2)==0
            str=['运行中...',num2str(kk/1e6),'秒'];
            waitbar(kk/20e6,h,str)
        end
        %% detect dual
        for ii=1:kk % 卷积
            xx(kk)=xx(kk)+(r_x(ii)+f_xq(ii))*g_xy(kk+1-ii);
            yy(kk)=yy(kk)+(r_y(ii)+f_yq(ii))*g_xy(kk+1-ii);
        end
        %% 解调.
        if mod((kk-1),8)==0
            kkT2=(kk-1)/8+1;
            for ii=1:kkT2 % 卷积
                lpf_amp(kkT2)=lpf_amp(kkT2)+(xx((ii-1)*8+1)*...
                    cos((delta_ome_r(idivide(kkT2-1,int32(125),'fix')+1)+...
                    ome0)*t((ii-1)*8+1)))*g_LPF1(kkT2+1-ii);
                lpf_LOOP(kkT2)=lpf_LOOP(kkT2)+(xx((ii-1)*8+1)*...
                    sin((delta_ome_r(idivide(kkT2-1,int32(125),'fix')+1)+...
                    ome0)*t((ii-1)*8+1)))*g_LPF1(kkT2+1-ii);
                lpf_yc(kkT2)=lpf_yc(kkT2)+(yy((ii-1)*8+1)*...
                    cos((delta_ome_r(idivide(kkT2-1,int32(125),'fix')+1)+...
                    ome0)*t((ii-1)*8+1)))*g_LPF1(kkT2+1-ii);
                lpf_ys(kkT2)=lpf_ys(kkT2)+(yy((ii-1)*8+1)*...
                    sin((delta_ome_r(idivide(kkT2-1,int32(125),'fix')+1)+...
                    ome0)*t((ii-1)*8+1)))*g_LPF1(kkT2+1-ii);
            end
        end
        %% 多PI/频率跟踪算法
        if mod((kk-1),1000)==0
            kkT3=(kk-1)/1000+1;
            PD(kkT3)=jdjs(lpf_LOOP(kkT2),lpf_LOOP(kkT2));
            for ii=1:kkT3-1 % 卷积
                f_xc(kkT3)=f_xc(kkT3)+(A0-lpf_amp((ii-1)*125+1))*(I1*1e-3);
                f_yc(kkT3)=f_yc(kkT3)+(-lpf_yc((ii-1)*125+1))*(I2*1e-3);
                f_ys(kkT3)=f_ys(kkT3)+(-lpf_ys((ii-1)*125+1))*(I3*1e-3);
                delta_ome_r(kkT3)=delta_ome_r(kkT3)+(-PD(ii))*(I3*1e-3);
            end
            f_xc(kkT3)=f_xc(kkT3)+(A0-lpf_amp((kkT3-1)*125+1))*(P1);
            f_yc(kkT3)=f_yc(kkT3)+(-lpf_yc((kkT3-1)*125+1))*(P2);
            f_ys(kkT3)=f_ys(kkT3)+(-lpf_ys((kkT3-1)*125+1))*(P3);
            delta_ome_r(kkT3)=delta_ome_r(kkT3)+(-PD(kkT3))*(P3);
        end
        %% 调制
        f_xq(kk)=f_xc(idivide(kk-1,int32(1000),'fix')+1)*...
            cos((delta_ome_r(idivide(kk-1,int32(1000),'fix')+1)+ome0)*t(kk)+pi/2);
        f_yq(kk)=f_yc(idivide(kk-1,int32(1000),'fix')+1)*...
            cos((delta_ome_r(idivide(kk-1,int32(1000),'fix')+1)+ome0)*t(kk)+pi/2)+...
            f_ys(idivide(kk-1,int32(1000),'fix')+1)*...
            sin((delta_ome_r(idivide(kk-1,int32(1000),'fix')+1)+ome0)*t(kk)+pi/2);        
        %% drive dual
        % drive voltage trans
    end
end
delete(h);
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