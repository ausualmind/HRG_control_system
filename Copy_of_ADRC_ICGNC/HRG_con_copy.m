function O = HRG_con(PIs,ADRC,preset,GHox,con_style) % HRG控制电路
O=zeros(4,1);
hh=waitbar(0,'please wait');
P1=PIs(1,1); I1=PIs(1,2); P2=PIs(2,1); I2=PIs(2,2);
P3=PIs(3,1); I3=PIs(3,2); P4=PIs(4,1); I4=PIs(4,2);
TT=100e6; xxx=zeros(4,TT+1); FF=zeros(2,TT+1); phi=zeros(1,TT+1);%0.01秒
x=zeros(1,TT+1); y=zeros(1,TT+1); p=zeros(1,TT+1); q=zeros(1,TT+1);
x(1)=10e-6*GHox.ome0; F1=zeros(1,TT+1); F2=zeros(1,TT+1);
pVc=zeros(1,TT/8+1); qVc=zeros(1,TT/8+1); % 解调值初始化
pVs=zeros(1,TT/8+1); qVs=zeros(1,TT/8+1);
Cp=zeros(1,TT/8+1); Cq=zeros(1,TT/8+1); % 滤波值初始化
Sp=zeros(1,TT/8+1); Sq=zeros(1,TT/8+1);
E=zeros(1,TT/1e3+1); % p轴振幅
S_delphi=zeros(1,TT/1e3+1); C_delphi=zeros(1,TT/1e3+1); % 组合运算频差正余弦值初始化
fome=zeros(1,TT/1e3+1); fps=zeros(1,TT/1e3+1); % 控制力初始化
fqc=zeros(1,TT/1e3+1); fqs=zeros(1,TT/1e3+1);
err_ome=zeros(1,TT/1e3+1); err_ps=zeros(1,TT/1e3+1); % 误差描述初始化
err_qc=zeros(1,TT/1e3+1); err_qs=zeros(1,TT/1e3+1);
xd12_ps=zeros(2,TT/1e3+1); z123_ps=zeros(3,TT/1e3+1); %ADRC所用状态量
xd12_qc=zeros(2,TT/1e3+1); z123_qc=zeros(3,TT/1e3+1);
xd12_qs=zeros(2,TT/1e3+1); z123_qs=zeros(3,TT/1e3+1);
for k=1:TT % 仿真步长
    kk=idivide(k-1,int32(8))+1; % 当k=1时 kk=1
    kkk=idivide(k-1,int32(1e3))+1; % 当k=1时 kkk=1
    if mod((k-1),1e5)==0
        str=['运行中...',num2str((k-1)/1e6),'秒'];
        waitbar(k/TT,hh,str)
    end
    xxx(:,k)=[x(k) y(k) p(k) q(k)]'; FF(:,k)=[F1(k) F2(k)]';
    x(k+1)=GHox.GG(1,:)*xxx(:,k) + GHox.HH(1,:)*FF(:,k);
    y(k+1)=GHox.GG(2,:)*xxx(:,k) + GHox.HH(2,:)*FF(:,k);
    p(k+1)=GHox.GG(3,:)*xxx(:,k) + GHox.HH(3,:)*FF(:,k);
    q(k+1)=GHox.GG(4,:)*xxx(:,k) + GHox.HH(4,:)*FF(:,k);
    %% 解调滤波
    if mod(k-1,8)==0
        pVc(kk)=2e5*p(k)*2*cos(phi(k)); qVc(kk)=2e5*q(k)*2*cos(phi(k));%正交
        pVs(kk)=2e5*p(k)*2*sin(phi(k)); qVs(kk)=2e5*q(k)*2*sin(phi(k));%同相
        %% 椭圆滤波器 采样频率125000Hz 通带频率150Hz
%         a0=1; % 滤波器分子参数
%         a1=-1.509819631626288405357172450749203562737;
%         a2=1;
%         b1=-1.991695518332145642403929741703905165195; % 滤波器分母参数
%         b2= 0.991757939352337358762667918199440464377;
%         GAIN=0.000113494534654183903119896437150515567;
        %% 椭圆滤波器 采样频率125000Hz 通带频率50Hz
        a0=1; % 滤波器分子参数
        a1=-1.938878010963012377132486108166631311178;
        a2=1;
        b1=-1.997238096155844466750295396195724606514; % 滤波器分母参数
        b2= 0.997245050970188895433921061339788138866;
        GAIN=0.000101411699891197794711275581569509541 ;
        %% 
        if kk>2 
        Cp(kk)=-b1*Cp(kk-1)-b2*Cp(kk-2)+GAIN*( a0*pVc(kk)+a1*pVc(kk-1)+a2*pVc(kk-2) );
        Cq(kk)=-b1*Cq(kk-1)-b2*Cq(kk-2)+GAIN*( a0*qVc(kk)+a1*qVc(kk-1)+a2*qVc(kk-2) );
        Sp(kk)=-b1*Sp(kk-1)-b2*Sp(kk-2)+GAIN*( a0*pVs(kk)+a1*pVs(kk-1)+a2*pVs(kk-2) );
        Sq(kk)=-b1*Sq(kk-1)-b2*Sq(kk-2)+GAIN*( a0*qVs(kk)+a1*qVs(kk-1)+a2*qVs(kk-2) );
        end
    end
    %% 多PI/多ADRC/频率跟踪算法
    if mod(k-1,1e3)==0
        Ts3=1/1e3;
        E(kkk)=sqrt(Cp(kk)^2+Sp(kk)^2)*1.122;
        S_delphi(kkk)=2*( Cp(kk)*Sp(kk)+Cq(kk)*Sq(kk) );
        C_delphi(kkk)=Cp(kk)^2+Cq(kk)^2-Sp(kk)^2-Sq(kk)^2;
        err_ome(kkk)=-1/2*jdjs(-S_delphi(kkk),-C_delphi(kkk)); % 误差评定指标！！
        if kkk>1
            fome(kkk)=fome(kkk-1)+P1*err_ome(kkk)+...
                (I1*Ts3-P1)*err_ome(kkk-1); % PLL中LP输出
        end
        err_ps(kkk)=preset.ps.xd-E(kkk); % 误差评定指标！！
        err_qc(kkk)=preset.qc.xd-Cq(kk)*1.122; % 误差评定指标！！
        err_qs(kkk)=preset.qs.xd-Sq(kk)*1.122; % 误差评定指标！！
        %% PI算法/ADRC算法选择
        if (kkk>1) && (con_style==1)
            fps(kkk)=fps(kkk-1)+P2*err_ps(kkk)+(I2*Ts3-P2)*err_ps(kkk-1);
%             fqc(kkk)=fqc(kkk-1)+P3*err_qc(kkk)+(I3*Ts3-P3)*err_qc(kkk-1);
%             fqs(kkk)=fqs(kkk-1)+P4*err_qs(kkk)+(I4*Ts3-P4)*err_qs(kkk-1);
        end
        %% ADRC算法起效
        if (kkk>1) && (con_style==2)
            [fps(kkk-1:kkk),xd12_ps(:,kkk-1:kkk),z123_ps(:,kkk-1:kkk)]=...
                ADRC_w(ADRC,Sp(kk),preset.ps,xd12_ps(:,kkk-1:kkk),z123_ps(:,kkk-1:kkk),...
                fps(kkk-1:kkk));
%             [fqc(kkk-1:kkk),xd12_qc(:,kkk-1:kkk),z123_qc(:,kkk-1:kkk)]=...
%                 ADRC_w(ADRC,Cq(kk),preset.qc,xd12_qc(:,kkk-1:kkk),z123_qc(:,kkk-1:kkk),...
%                 fqc(kkk-1:kkk));
%             [fqs(kkk-1:kkk),xd12_qs(:,kkk-1:kkk),z123_qs(:,kkk-1:kkk)]=...
%                 ADRC_w(ADRC,Sq(kk),preset.qs,xd12_qs(:,kkk-1:kkk),z123_qs(:,kkk-1:kkk),...
%                 fqs(kkk-1:kkk));
        end
    end
    %%% 调制
    phi(k+1)=phi(k)+(fome(kkk)+GHox.ome0)/1e6;
    F1(k+1)=1/GHox.xx* fps(kkk)* sin(phi(k+1)+pi/2);
%     F2(k+1)=1/GHox.xx* ( fqs(kkk)* sin(phi(k+1)+pi/2)+fqc(kkk)* cos(phi(k+1)+pi/2) );
end
%% 计算误差性能指标
for iii=1:kkk
    O(1)=O(1)+double(iii-1)*Ts3^2*abs( err_ome(iii) );
    O(2)=O(2)+double(iii-1)*Ts3^2*abs( err_ps(iii) );
    O(3)=O(3)+double(iii-1)*Ts3^2*abs( err_qc(iii) );
    O(4)=O(4)+double(iii-1)*Ts3^2*abs( err_qs(iii) );
end
delete(hh);
end
function y = jdjs(S,C)
    y = 0; 
    if S >= 0
        if C > 0, y = atan(S/C); end
        if C < 0, y = (pi + atan(S/C) ); end
        if C == 0, y = pi/2; end
    else
        if C > 0, y = ( atan(S/C) ); end
        if C < 0, y = ( atan(S/C) - pi); end
        if C == 0, y = -pi/2; end
    end
end
