function [Swarm,cg_curve] = PSO_w(inref,fobj,PIs,ADRC,preset,GHox,con_style) % HRG控制电路
ite_con=zeros(inref.Max_iteration,inref.noP); % 杂交变异判定数组/并判断是否满足终止条件
err_select1234=2; % 1表示J1/ 2表示J2/ 3表示J3/ 4表示J4/
%% 初始化赋值（1变量内存配置，2等高线绘制，3边界描述，4粒子分布，5初始粒子分布图）
X1max=inref.X1max; X1min=inref.X1min; % X1/X2参数上下界
X2max=inref.X2max; X2min=inref.X2min;
noP=inref.noP; nn=inref.nn;
% PSO paramters
Max_iteration = inref.Max_iteration;
Vmax=inref.Vmax;
wmax=inref.wmax; wmin=inref.wmin;
c1=inref.c1; c2=inref.c2;
%% 1变量内存配置
first_vel=zeros(inref.noP,2); second_vel=zeros(inref.noP,2);
first_loc=zeros(inref.noP,2); second_loc=zeros(inref.noP,2);
moving_x=zeros(inref.noP,inref.nn); moving_y=zeros(inref.noP,inref.nn); % 行数：粒子数，列数：渐进区间
moving_vx=zeros(inref.noP,inref.nn); moving_vy=zeros(inref.noP,inref.nn);
h2=zeros(inref.noP);
%% 2等高线绘制
figure
% Renderer用于屏幕和图片的渲染模式/有效值：painters、zbuffer、OpenGL/缺省值：系统自动选择
set(gcf,'Renderer','OpenGL');
x = linspace(X1min,X1max,5);
y = linspace(X2min,X2max,5);
[x_new, y_new] = meshgrid(x,y);
z=zeros(size(x_new, 1),size(x_new , 2));
hh=waitbar(0,'等高线绘制');
for k1 = 1: size(x_new, 1)
    for k2 = 1 : size(x_new , 2)
        str=['运行中...',num2str((k1-1)*5+k2),'个点'];
        waitbar(((k1-1)*5+k2)/25,hh,str)
        X = [ x_new(k1,k2) , y_new(k1, k2) ];
%         PIs(GHox.err1234,:)=X; % PI参数寻优
        if err_select1234==1
%             ADRC.beta01=X(1); ADRC.beta02=X(2); % 对ADRC分离性原则求参
%             ADRC.beta01=X(1); ADRC.beta03=X(2);
            ADRC.beta02=X(1); ADRC.beta03=X(2);
        end
        if (err_select1234==2)||(err_select1234==3)
            ADRC.beta1=X(1); ADRC.beta2=X(2);
        end
        J = fobj( PIs,ADRC,preset,GHox,con_style );
        z(k1,k2) =J(err_select1234);
    end
end
h = contour(x_new, y_new, z , 5000);
set(gca,'XLimMode','manual');
set(gca,'YLimMode','manual');
axis([X1min,X1max,X2min,X2max])
hold on
view(2)
shading interp
%set(h,'EraseMode','normal')
%% 3边界描述
nVar = 2; % 每个粒子的坐标维度
% Objective function details
lb = [X1min,X2min];
ub = [X1max,X2max];

cg_curve = zeros(1,inref.Max_iteration);

x = X1min:(X1max-X1min)/5:X1max;
y = X2min:(X2max-X2min)/5:X2max;
%% 4粒子分布
idx = 1;
show_vel = 1;
for t1 = 1: length (x) % inref.noP个粒子的初始位置
    for t2 =  1: length (y)
        Swarm.Particles(idx).X = [x(t1) y(t2)];
        idx = idx+1;
    end
end

%% 5初始粒子分布图
hh=waitbar(0,'初始化 粒子数');
for k = 1: inref.noP
    str=['运行中...',num2str(k),'个'];
    waitbar(k/inref.noP,hh,str)
    
    Swarm.Particles(k).V = rand(1,nVar); % 随机数
    Swarm.Particles(k).PBEST.X = rand(1,nVar).*inref.Vmax; % 随机数
    Swarm.Particles(k).PBEST.O= inf; % for minimization problems
    h(k) = plot(Swarm.Particles(k).X(1),Swarm.Particles(k).X(2),'ok', 'markerFaceColor','k');
%     set(h(k),'EraseMode','normal')
    if show_vel == 1
        p1 = [Swarm.Particles(k).X(1) Swarm.Particles(k).X(2)]; % First Point
        p2 = [Swarm.Particles(k).V(1) Swarm.Particles(k).V(2)]; % Second Point
        dp = p2-p1/100;                         % Difference
        
        h2(k) = plot([p1(1) dp(1)],[p1(2) dp(2)],'-k' , 'LineWidth',2.5);
%         set(h2(k),'EraseMode','normal')
    end
end
delete(hh);

Swarm.GBEST.X=zeros(1,nVar);
Swarm.GBEST.O= inf;
%% 开始迭代
hh=waitbar(0,'迭代次数');
for t=1:inref.Max_iteration % main loop
    str=['运行中...',num2str(t),'次'];
    waitbar(t/inref.Max_iteration,hh,str)
    %% 更新粒子历史最优和全局最优
    for k=1:inref.noP % Calculate objective function for each particle
%         PIs(GHox.err1234,:)=Swarm.Particles(k).X;
        % 对ADRC分离性原则求参
        if err_select1234==1
%         ADRC.beta01=Swarm.Particles(k).X(1); ADRC.beta02=Swarm.Particles(k).X(2);
%         ADRC.beta01=Swarm.Particles(k).X(1); ADRC.beta03=Swarm.Particles(k).X(2);
            ADRC.beta02=Swarm.Particles(k).X(1); ADRC.beta03=Swarm.Particles(k).X(2);
        end
        if (err_select1234==2)||(err_select1234==3)        
            ADRC.beta1=Swarm.Particles(k).X(1); ADRC.beta2=Swarm.Particles(k).X(2);
        end
        J=fobj( PIs,ADRC,preset,GHox,con_style );
        Swarm.Particles(k).O=J(err_select1234);
        
        if(Swarm.Particles(k).O < Swarm.Particles(k).PBEST.O) % 粒子本身经历过的最优位置
            Swarm.Particles(k).PBEST.O = Swarm.Particles(k).O;
            Swarm.Particles(k).PBEST.X = Swarm.Particles(k).X;
        end
        if(Swarm.Particles(k).O < Swarm.GBEST.O) % 粒子群经历过的最优位置
            Swarm.GBEST.O = Swarm.Particles(k).O;
            Swarm.GBEST.X = Swarm.Particles(k).X;
            % k个粒子的更新数组，若当前迭代中存在粒子最优值更新，则数组中一些元素变为0，否则全为1
            ite_con(t,k)=0;
        else
            ite_con(t,k)=1;
        end
    end
    %% PSO杂交变异
    if (t>=5)&&(sum(ite_con(t-4:t,:),'all')==5*36)
        select=round(rand(1,20)*35)+1;
        for s_i=1:20
            r=rand(1,nVar);
            Swarm.Particles(select(s_i)).X=r.*Swarm.Particles(select(s_i)).PBEST.X+...
                (1-r).*Swarm.Particles(select(s_i)).X;
        end
    end
    %% Update the inertia weight
    w=inref.wmax-t*((inref.wmax-inref.wmin)/inref.Max_iteration); % 更新惯性权重w
    %% Update the Velocity and Position of particles
    for k=1:inref.noP % 适应度函数值计算
        first_vel(k,:) = Swarm.Particles(k).V;
        Swarm.Particles(k).V  = w .* Swarm.Particles(k).V + ...  % inertia
            inref.c1 .* rand(1,nVar) .* (Swarm.Particles(k).PBEST.X - Swarm.Particles(k).X ) +  ...   % congnitive
            inref.c2 .* rand(1,nVar).* (Swarm.GBEST.X - Swarm.Particles(k).X) ;  % social
        
        second_vel(k,:) = Swarm.Particles(k).V;
        
        index = find(Swarm.Particles(k).V > inref.Vmax);
        Swarm.Particles(k).V(index) = inref.Vmax(index)*rand;
        
        index = find(Swarm.Particles(k).V < -inref.Vmax);
        Swarm.Particles(k).V(index) = -inref.Vmax(index)*rand;
        
        first_loc(k,:) = Swarm.Particles(k).X ;
        Swarm.Particles(k).X = Swarm.Particles(k).X + Swarm.Particles(k).V;
        second_loc(k,:) = Swarm.Particles(k).X ;
        
        indx = find( Swarm.Particles(k).X  > ub);
        Swarm.Particles(k).X (indx) = ub(indx);
        
        indx = find( Swarm.Particles(k).X  < lb);
        Swarm.Particles(k).X (indx) = lb(indx);
        
        moving_x(k,:) =  linspace(first_loc(k,1),second_loc(k,1),inref.nn);
        moving_y(k,:) = linspace(first_loc(k,2),second_loc(k,2),inref.nn);
        
        moving_vx(k,:) =  linspace(first_vel(k,1),second_vel(k,1),inref.nn);
        moving_vy(k,:) = linspace(first_vel(k,2),second_vel(k,2),inref.nn);
    end
    %% 动画演示每个粒子从上一位置到目前位置的运动曲线
    for rr = 1: inref.nn % inref.nn个渐进区间
        for JJ = 1:inref.noP % inref.noP个粒子
            set(h(JJ),'XData', moving_x(JJ,rr),'YData', moving_y(JJ,rr));
            if show_vel == 1
                p1 = [moving_x(JJ,rr) moving_y(JJ,rr)]; % First Point
                p2 = [moving_x(JJ,end) moving_y(JJ,end)];
                dp = (p2+p1)/2;                         % Difference
                dp = (p1+dp)/2;
                set(h2(JJ),'XData', [p1(1), dp(1)],'YData',  [p1(2), dp(2)]);
            end
        end
        drawnow limitrate nocallbacks
    end
    cg_curve(t) = Swarm.GBEST.O;
    %% 判断是否满足迭代终止条件
    if (t>=10)&&(sum(ite_con(t-9:t,:),'all')==10*36)
        break;
    end
end
%% 绘制fobj函数值随迭代次数变化图
delete(hh);
figure
plot(cg_curve)
xlabel('Iteration')
