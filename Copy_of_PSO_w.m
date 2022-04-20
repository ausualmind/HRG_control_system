function [Swarm,cg_curve] = PSO_w(inref,fobj,PIs,ADRC,preset,GHox,con_style) % HRG���Ƶ�·
ite_con=zeros(inref.Max_iteration,inref.noP); % �ӽ������ж�����/���ж��Ƿ�������ֹ����
err_select1234=2; % 1��ʾJ1/ 2��ʾJ2/ 3��ʾJ3/ 4��ʾJ4/
%% ��ʼ����ֵ��1�����ڴ����ã�2�ȸ��߻��ƣ�3�߽�������4���ӷֲ���5��ʼ���ӷֲ�ͼ��
X1max=inref.X1max; X1min=inref.X1min; % X1/X2�������½�
X2max=inref.X2max; X2min=inref.X2min;
noP=inref.noP; nn=inref.nn;
% PSO paramters
Max_iteration = inref.Max_iteration;
Vmax=inref.Vmax;
wmax=inref.wmax; wmin=inref.wmin;
c1=inref.c1; c2=inref.c2;
%% 1�����ڴ�����
first_vel=zeros(inref.noP,2); second_vel=zeros(inref.noP,2);
first_loc=zeros(inref.noP,2); second_loc=zeros(inref.noP,2);
moving_x=zeros(inref.noP,inref.nn); moving_y=zeros(inref.noP,inref.nn); % ����������������������������
moving_vx=zeros(inref.noP,inref.nn); moving_vy=zeros(inref.noP,inref.nn);
h2=zeros(inref.noP);
%% 2�ȸ��߻���
figure
% Renderer������Ļ��ͼƬ����Ⱦģʽ/��Чֵ��painters��zbuffer��OpenGL/ȱʡֵ��ϵͳ�Զ�ѡ��
set(gcf,'Renderer','OpenGL');
x = linspace(X1min,X1max,5);
y = linspace(X2min,X2max,5);
[x_new, y_new] = meshgrid(x,y);
z=zeros(size(x_new, 1),size(x_new , 2));
hh=waitbar(0,'�ȸ��߻���');
for k1 = 1: size(x_new, 1)
    for k2 = 1 : size(x_new , 2)
        str=['������...',num2str((k1-1)*5+k2),'����'];
        waitbar(((k1-1)*5+k2)/25,hh,str)
        X = [ x_new(k1,k2) , y_new(k1, k2) ];
%         PIs(GHox.err1234,:)=X; % PI����Ѱ��
        if err_select1234==1
%             ADRC.beta01=X(1); ADRC.beta02=X(2); % ��ADRC������ԭ�����
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
%% 3�߽�����
nVar = 2; % ÿ�����ӵ�����ά��
% Objective function details
lb = [X1min,X2min];
ub = [X1max,X2max];

cg_curve = zeros(1,inref.Max_iteration);

x = X1min:(X1max-X1min)/5:X1max;
y = X2min:(X2max-X2min)/5:X2max;
%% 4���ӷֲ�
idx = 1;
show_vel = 1;
for t1 = 1: length (x) % inref.noP�����ӵĳ�ʼλ��
    for t2 =  1: length (y)
        Swarm.Particles(idx).X = [x(t1) y(t2)];
        idx = idx+1;
    end
end

%% 5��ʼ���ӷֲ�ͼ
hh=waitbar(0,'��ʼ�� ������');
for k = 1: inref.noP
    str=['������...',num2str(k),'��'];
    waitbar(k/inref.noP,hh,str)
    
    Swarm.Particles(k).V = rand(1,nVar); % �����
    Swarm.Particles(k).PBEST.X = rand(1,nVar).*inref.Vmax; % �����
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
%% ��ʼ����
hh=waitbar(0,'��������');
for t=1:inref.Max_iteration % main loop
    str=['������...',num2str(t),'��'];
    waitbar(t/inref.Max_iteration,hh,str)
    %% ����������ʷ���ź�ȫ������
    for k=1:inref.noP % Calculate objective function for each particle
%         PIs(GHox.err1234,:)=Swarm.Particles(k).X;
        % ��ADRC������ԭ�����
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
        
        if(Swarm.Particles(k).O < Swarm.Particles(k).PBEST.O) % ���ӱ�������������λ��
            Swarm.Particles(k).PBEST.O = Swarm.Particles(k).O;
            Swarm.Particles(k).PBEST.X = Swarm.Particles(k).X;
        end
        if(Swarm.Particles(k).O < Swarm.GBEST.O) % ����Ⱥ������������λ��
            Swarm.GBEST.O = Swarm.Particles(k).O;
            Swarm.GBEST.X = Swarm.Particles(k).X;
            % k�����ӵĸ������飬����ǰ�����д�����������ֵ���£���������һЩԪ�ر�Ϊ0������ȫΪ1
            ite_con(t,k)=0;
        else
            ite_con(t,k)=1;
        end
    end
    %% PSO�ӽ�����
    if (t>=5)&&(sum(ite_con(t-4:t,:),'all')==5*36)
        select=round(rand(1,20)*35)+1;
        for s_i=1:20
            r=rand(1,nVar);
            Swarm.Particles(select(s_i)).X=r.*Swarm.Particles(select(s_i)).PBEST.X+...
                (1-r).*Swarm.Particles(select(s_i)).X;
        end
    end
    %% Update the inertia weight
    w=inref.wmax-t*((inref.wmax-inref.wmin)/inref.Max_iteration); % ���¹���Ȩ��w
    %% Update the Velocity and Position of particles
    for k=1:inref.noP % ��Ӧ�Ⱥ���ֵ����
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
    %% ������ʾÿ�����Ӵ���һλ�õ�Ŀǰλ�õ��˶�����
    for rr = 1: inref.nn % inref.nn����������
        for JJ = 1:inref.noP % inref.noP������
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
    %% �ж��Ƿ����������ֹ����
    if (t>=10)&&(sum(ite_con(t-9:t,:),'all')==10*36)
        break;
    end
end
%% ����fobj����ֵ����������仯ͼ
delete(hh);
figure
plot(cg_curve)
xlabel('Iteration')
