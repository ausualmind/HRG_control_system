function J = HRG_con(PIs,ADRC,preset,GHox,con_style) % HRG控制电路
O=zeros(4,1); J=zeros(4,1);
% hh=waitbar(0,'please wait');
P1=PIs(1,1); I1=PIs(1,2); P2=PIs(2,1); I2=PIs(2,2);
P3=PIs(3,1); I3=PIs(3,2); P4=PIs(4,1); I4=PIs(4,2);
TT=10*40e3; beishu=40; ts1=40e3; % 0.01秒
xxx=zeros(4,TT+1); FF=zeros(2,TT+1); phi=zeros(1,TT+1);
x=zeros(1,TT+1); y=zeros(1,TT+1); p=zeros(1,TT+1); q=zeros(1,TT+1);
x(1)=10e-6*GHox.ome0; F1=zeros(1,TT+1); F2=zeros(1,TT+1);
pVc=zeros(1,TT+1); qVc=zeros(1,TT+1); % 解调值初始化
pVs=zeros(1,TT+1); qVs=zeros(1,TT+1);
a=zeros(1,TT+1); b=zeros(1,TT+1); % 滤波值初始化
m=zeros(1,TT+1); n=zeros(1,TT+1);
a_real=zeros(1,TT/beishu+1); b_real=zeros(1,TT/beishu+1);
m_real=zeros(1,TT/beishu+1); n_real=zeros(1,TT/beishu+1);
E=zeros(1,TT/beishu+1); % p轴振幅
S_delphi=zeros(1,TT/beishu+1); C_delphi=zeros(1,TT/beishu+1);
precession=zeros(1,TT/beishu+1); pre=zeros(1,TT/beishu+1); % 进动角解算
pre_show=zeros(1,TT/beishu+1);
fome=zeros(1,TT/beishu+1); fps=zeros(1,TT/beishu+1); % 控制力初始化
fps1=zeros(1,TT/beishu+1);
fqc=zeros(1,TT/beishu+1); fqs=zeros(1,TT/beishu+1);
err_ome=zeros(1,TT/beishu+1); err_ps=zeros(1,TT/beishu+1); % 误差描述初始化
err_qc=zeros(1,TT/beishu+1); err_qs=zeros(1,TT/beishu+1);
xd12_ps=zeros(2,TT/beishu+1); z123_ps=zeros(3,TT/beishu+1); %ADRC所用状态量
xd12_qc=zeros(2,TT/beishu+1); z123_qc=zeros(3,TT/beishu+1);
xd12_qs=zeros(2,TT/beishu+1); z123_qs=zeros(3,TT/beishu+1);
for k=1:TT+1 % 仿真步长
%     kk=idivide(k-1,int32(8))+1; % 当k=1时 kk=1
    kk=k;
%     kkk=idivide(k-1,int32(1e3))+1; % 当k=1时 kkk=1
    kkk=idivide(k-1,int32(beishu))+1; % 当k=1时 kkk=1
%     if mod((k-1),ts1)==0
%         str=['运行中...',num2str((k-1)/ts1),'秒'];
%         waitbar(k/TT,hh,str)
%     end
    xxx(:,k)=[x(k) y(k) p(k) q(k)]'; FF(:,k)=[F1(k) F2(k)]';
    x(k+1)=GHox.GG(1,:)*xxx(:,k) + GHox.HH(1,:)*FF(:,k);
    y(k+1)=GHox.GG(2,:)*xxx(:,k) + GHox.HH(2,:)*FF(:,k);
    p(k+1)=GHox.GG(3,:)*xxx(:,k) + GHox.HH(3,:)*FF(:,k);
    q(k+1)=GHox.GG(4,:)*xxx(:,k) + GHox.HH(4,:)*FF(:,k);
    %% 解调滤波
%     if mod(k-1,20)==0
        pVc(kk)=2e5*p(k)*2*cos(phi(k)); qVc(kk)=2e5*q(k)*2*cos(phi(k));%正交
        pVs(kk)=2e5*p(k)*2*sin(phi(k)); qVs(kk)=2e5*q(k)*2*sin(phi(k));%同相
        %% 椭圆滤波器 采样频率4万Hz 通带频率50Hz
        a0=1; % 滤波器分子参数
        a1=-1.473607447387377034075939263857435435057;
        a2=1;
        b1=-1.991348282585754470019878681341651827097 ; % 滤波器分母参数
        b2= 0.991416002063946311118058929423568770289;
        GAIN=0.000114657869244153337779691670395010306 ;
        %% 
        if kk>2
        a(kk)=-b1*a(kk-1)-b2*a(kk-2)+GAIN*( a0*pVc(kk)+a1*pVc(kk-1)+a2*pVc(kk-2) );
        b(kk)=-b1*b(kk-1)-b2*b(kk-2)+GAIN*( a0*qVc(kk)+a1*qVc(kk-1)+a2*qVc(kk-2) );
        m(kk)=-b1*m(kk-1)-b2*m(kk-2)+GAIN*( a0*pVs(kk)+a1*pVs(kk-1)+a2*pVs(kk-2) );
        n(kk)=-b1*n(kk-1)-b2*n(kk-2)+GAIN*( a0*qVs(kk)+a1*qVs(kk-1)+a2*qVs(kk-2) );
        end
%     end
    %% 多PI/多ADRC/频率跟踪算法
    if mod(k-1,beishu)==0
        Ts3=1/1e3;
        a_real(kkk)=a(kk)*1.122/2;
        b_real(kkk)=b(kk)*1.122/2;
        m_real(kkk)=m(kk)*1.122/2;
        n_real(kkk)=n(kk)*1.122/2;
        E(kkk)=sqrt( a(kk)^2+b(kk)^2+m(kk)^2+n(kk)^2 )*1.122;
        S_delphi(kkk)=2*( a(kk)*m(kk)+b(kk)*n(kk) );
        C_delphi(kkk)=a(kk)^2+b(kk)^2-m(kk)^2-n(kk)^2;
        err_ome(kkk)=-1/2*pre_solve(-S_delphi(kkk),-C_delphi(kkk)); % 误差评定指标！！
        if kkk>1
            fome(kkk)=fome(kkk-1)+P1*err_ome(kkk)+...
                (I1*Ts3-P1)*err_ome(kkk-1); % PLL中LP输出
            precession(kkk)=1/4*pre_solve( 2*( a(kk)*b(kk)+m(kk)*n(kk) ),...
                a(kk)^2+m(kk)^2-b(kk)^2-n(kk)^2 );
            e1=( precession(kkk)-precession(kkk-1) )/pi*180;
            e2=precession(kkk)/pi*180-pre(kkk-1);
            if abs(e1)>80
                if ( abs(e2)+10 )<abs(e1)
                    pre(kkk)=precession(kkk)/pi*180;
                else
                    pre(kkk)=precession(kkk)/pi*180-idivide(e2,int32(90),'round')*90;
                end
            else
                pre(kkk)=pre(kkk-1)+ e1;
            end
        end
        pre_show(kkk)=pre_show_v(pre(kkk)/180*pi); %进动角展示
        err_ps(kkk)=( preset.ps.xd-E(kkk) ); % 误差评定指标！！
        err_qc(kkk)=preset.qc.xd-b(kk)*1.122; % 误差评定指标！！
        err_qs(kkk)=preset.qs.xd-n(kk)*1.122; % 误差评定指标！！
        %% PI算法/ADRC算法选择
        if (kkk>1) && (con_style==1)
            fps(kkk)=fps(kkk-1)+P2*err_ps(kkk)+...
                (I2*Ts3-P2)*( preset.ps.xd-E(kkk-1) );
%             fqc(kkk)=fqc(kkk-1)+P3*err_qc(kkk)+(I3*Ts3-P3)*err_qc(kkk-1);
%             fqs(kkk)=fqs(kkk-1)+P4*err_qs(kkk)+(I4*Ts3-P4)*err_qs(kkk-1);
        end
        %% ADRC算法起效
        if (kkk>1) && (con_style==2)
            [fps(kkk-1:kkk),xd12_ps(:,kkk-1:kkk),z123_ps(:,kkk-1:kkk)]=...
                ADRC_w(ADRC,E(kkk),preset.ps,xd12_ps(:,kkk-1:kkk),z123_ps(:,kkk-1:kkk),...
                fps(kkk-1:kkk));
%             [fqc(kkk-1:kkk),xd12_qc(:,kkk-1:kkk),z123_qc(:,kkk-1:kkk)]=...
%                 ADRC_w(ADRC,b(kk),preset.qc,xd12_qc(:,kkk-1:kkk),z123_qc(:,kkk-1:kkk),...
%                 fqc(kkk-1:kkk));
%             [fqs(kkk-1:kkk),xd12_qs(:,kkk-1:kkk),z123_qs(:,kkk-1:kkk)]=...
%                 ADRC_w(ADRC,n(kk),preset.qs,xd12_qs(:,kkk-1:kkk),z123_qs(:,kkk-1:kkk),...
%                 fqs(kkk-1:kkk));
        end
    end
    %%% 调制
    phi(k+1)=phi(k)+(fome(kkk)+GHox.ome0)/ts1;
    F1(k+1)=1/GHox.xx* fps(kkk)*sin(phi(k+1)+pi/2);
%     F2(k+1)=1/GHox.xx* ( fqs(kkk)* sin(phi(k+1)+pi/2)+fqc(kkk)* cos(phi(k+1)+pi/2) );
    if mod(k-1,ts1)==0
        %% 检查程序行
        aaaaaa_1=E(kkk);
        aaaaaa_2=err_ps(kkk);
        aaaaaa_3=fps(kkk);
        aaaaaa_4=1;
    end
end
%% 计算误差性能指标
for iii=1:kkk
    O(1)=O(1)+double(iii-1)*Ts3^2*abs( err_ome(iii) );
    O(2)=O(2)+double(iii-1)*Ts3^2*abs( err_ps(iii) );
    O(3)=O(3)+double(iii-1)*Ts3^2*abs( err_qc(iii) );
    O(4)=O(4)+double(iii-1)*Ts3^2*abs( err_qs(iii) );
end
for iii=1:kkk
    J(1)=J(1)+double(iii-1)*Ts3^2*abs( z123_ps(1,iii)-E(iii) );%ESO误差性能指标
    J(2)=J(2)+double(iii-1)*Ts3^2*abs( xd12_ps(1,iii)-z123_ps(1,iii) );%NLSEF误差性能指标1
    J(3)=J(3)+double(iii-1)*Ts3^2*abs( xd12_ps(2,iii)-z123_ps(2,iii) );%NLSEF误差性能指标2
end
J(4)=O(GHox.err1234);
% delete(hh);
end
function y = pre_solve(S,C)
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
function y = pre_show_v(u)
if sin(u) >= 0
    if cos(u) >= 0, y = asin(sin(u));
    else, y = pi - asin(sin(u)); end
else
    if cos(u) >= 0, y = ( asin(sin(u)) );
    else, y = - pi - asin(sin(u)); end
end
y = y/pi*180;
end
