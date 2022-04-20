clear all
close all
clc
tic
noP=36; nn=4;
first_vel=zeros(noP,2); second_vel=zeros(noP,2);
first_loc=zeros(noP,2); second_loc=zeros(noP,2);
moving_x=zeros(noP,nn); moving_y=zeros(noP,nn); % 行数：粒子数，列数：渐进区间
moving_vx=zeros(noP,nn); moving_vy=zeros(noP,nn);
h2=zeros(noP);
figure
% Renderer用于屏幕和图片的渲染模式/有效值：painters、zbuffer、OpenGL/缺省值：系统自动选择
set(gcf,'Renderer','OpenGL'); 
x = linspace(-10,10,30);
y = linspace(-10,10,30);
[x_new, y_new] = meshgrid(x,y);
z=zeros(size(x_new, 1),size(x_new , 2));
for k1 = 1: size(x_new, 1)
    for k2 = 1 : size(x_new , 2)
        X = [ x_new(k1,k2) , y_new(k1, k2) ];
        z(k1,k2) = ObjectiveFunction( X );
    end
end
h = contour(x_new, y_new, z , 20); % 等高线绘制
set(gca,'XLimMode','manual');
set(gca,'YLimMode','manual');
axis([-10,10,-10,10])
hold on
view(2)
shading interp
%set(h,'EraseMode','normal')
nVar = 2; % 每个粒子的坐标维度
% Objective function details
fobj = @ObjectiveFunction;
lb = -10 * ones(1,nVar) ;
ub = 10 * ones(1,nVar);
% PSO paramters
Max_iteration = 200;
Vmax=6;
wMax=0.9;
wMin=0.2;
c1=2;
c2=2;

cg_curve = zeros(1,Max_iteration);

x = -10:4:10;
y= x;

idx = 1;
show_vel = 1;
for t1 = 1: length (x) % noP个粒子的初始位置
    for t2 =  1: length (y)
        Swarm.Particles(idx).X = [x(t1) y(t2)];
        idx = idx+1;
    end
end

% Initializations
hh=waitbar(0,'初始化 粒子数');
for k = 1: noP
    str=['运行中...',num2str(k),'个'];
    waitbar(k/noP,hh,str)
    
    Swarm.Particles(k).V = rand(1,nVar); % 随机数
    Swarm.Particles(k).PBEST.X = rand(1,nVar)*Vmax; % 随机数
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

hh=waitbar(0,'迭代次数');
for t=1:Max_iteration % main loop
    str=['运行中...',num2str(t),'次'];
    waitbar(t/Max_iteration,hh,str)

    for k=1:noP % 更新粒子历史最优和全局最优
        %Calculate objective function for each particle
        Swarm.Particles(k).O=fobj( Swarm.Particles(k).X );
        
        if(Swarm.Particles(k).O < Swarm.Particles(k).PBEST.O) % 粒子本身经历过的最优位置
            Swarm.Particles(k).PBEST.O = Swarm.Particles(k).O;
            Swarm.Particles(k).PBEST.X = Swarm.Particles(k).X;
        end
        if(Swarm.Particles(k).O < Swarm.GBEST.O) % 粒子群经历过的最优位置
            Swarm.GBEST.O = Swarm.Particles(k).O;
            Swarm.GBEST.X = Swarm.Particles(k).X;
        end
    end
    %Update the inertia weight
    w=wMax-t*((wMax-wMin)/Max_iteration); % 更新惯性权重w
    %Update the Velocity and Position of particles
    for k=1:noP
        first_vel(k,:) = Swarm.Particles(k).V; % 适应度函数值计算
        Swarm.Particles(k).V  = w .* Swarm.Particles(k).V + ...  % inertia
            c1 .* rand(1,nVar) .* (Swarm.Particles(k).PBEST.X - Swarm.Particles(k).X ) +  ...   % congnitive
            c2 .* rand(1,nVar).* (Swarm.GBEST.X - Swarm.Particles(k).X) ;  % social
        
        second_vel(k,:) = Swarm.Particles(k).V;
        
        index = find(Swarm.Particles(k).V > Vmax);
        Swarm.Particles(k).V(index) = Vmax*rand;
        
        index = find(Swarm.Particles(k).V < -Vmax);
        Swarm.Particles(k).V(index) = -Vmax*rand;
        
        first_loc(k,:) = Swarm.Particles(k).X ;
        Swarm.Particles(k).X = Swarm.Particles(k).X + Swarm.Particles(k).V;
        second_loc(k,:) = Swarm.Particles(k).X ;
        
        indx = find( Swarm.Particles(k).X  > ub);
        Swarm.Particles(k).X (indx) = ub(1);
        
        indx = find( Swarm.Particles(k).X  < lb);
        Swarm.Particles(k).X (indx) = lb(1);
        
        moving_x(k,:) =  linspace(first_loc(k,1),second_loc(k,1),nn);
        moving_y(k,:) = linspace(first_loc(k,2),second_loc(k,2),nn);
        
        moving_vx(k,:) =  linspace(first_vel(k,1),second_vel(k,1),nn);
        moving_vy(k,:) = linspace(first_vel(k,2),second_vel(k,2),nn);
    end
    for rr = 1: nn % nn个渐进区间
        for JJ = 1:noP % noP个粒子
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
end
delete(hh);
figure
plot(cg_curve)
xlabel('Iteration')
toc