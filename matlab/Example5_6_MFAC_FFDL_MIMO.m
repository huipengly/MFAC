% clear all;
% close all;
clc
%%%%%%%%%%% FFDL based MFAC for TITO nonlinear system 2010-06-10 %%%%%%%%%%%%%
%采样时间
N=800;My=2;Mu=2;
%被控对象初始值
y=cell(1);dy=cell(1);u=cell(1);du=cell(1);
y{1}(1)=0;
y{1}(2)=0;
y{2}(1)=1;
y{2}(2)=1;
y{3}=y{1};
y1(1:3)=0,y1(2)=1,y2(1:3)=0,y2(2)=1; 

u{1}=zeros(Mu,1);u{2}=zeros(Mu,1);
%控制器阶数
Lu=3;
Ly=1;
eta=0.5;
miu=1;
rou=0.5;
lamda=1;
M=50;
epsilon=10^(-5);
%控制器伪偏导数初值
Fai=cell(1);Fai1=cell(1);

Fai{1}(1:My,1:Mu)=[0.1 0; 0 0.1];
Fai{1}(1:My,Mu+1:Mu*Lu)=0;
Fai{1}(1:My,Mu*Lu:Mu*Lu+My*Ly)=0;
Fai{2}=Fai{1};Fai{3}=Fai{1};

%期望值
yd=cell(1);
for k=1:N+1
    yd1(k)=5*sin(k/50)+2*cos(k/20);
    yd2(k)=2*sin(k/50)+5*cos(k/20);
    yd{k}(1)=yd1(k);
    yd{k}(2)=yd2(k);
end

%程序循环
for k=3:N
    %improved projection algorithm proposed by Hou (1999)
    zerou=zeros(Mu,1);zeroy=zeros(My,1);
    if k==2
        dH=u{k-1}(:)-0;
    else
        dH=u{k-1}(:)-u{k-2}(:);
    end
    for i=2:Lu
        if k>i+1
            dH=[dH;u{k-i}(:)-u{k-i-1}(:)];
        else
            dH=[dH;zerou];
        end
    end
    for i=1:Ly
        if k>i+1
            dH=[dH;y{k-i}(:)-y{k-i-1}(:)];
        else
            dH=[dH;zeroy];
        end
    end
    Fai{k}=Fai{k-1}+eta*(y{k}(:)-y{k-1}(:)-Fai{k-1}(:,:)*dH)*dH'/(miu+norm(dH,2)^2);
%%%%%%%%%% Reset algorithm %%%%%%%%%%%%%%%%
    Fai1{k}=Fai{k}(1:My,1:Mu);
    if (Fai1{k}(1,1)<epsilon)|(Fai1{k}(1,1)>M)
        Fai1{k}(1,1)=Fai1{3}(1,1);
    end
    if (Fai1{k}(2,2)<epsilon)|(Fai1{k}(2,2)>M)
        Fai1{k}(2,2)=Fai1{3}(2,2);
    end
    Fai{k}(1:My,1:Mu)=Fai1{k};
%   control law
    if k==2
        tmp=u{k-1}(:)-0;
    else
        tmp=u{k-1}(:)-u{k-2}(:);
    end
    for i=2:Lu-1
        if k>i+1
            tmp=[tmp;u{k-i}(:)-u{k-i-1}(:)];
        else
            tmp=[tmp;zerou];
        end    
    end
    for i=1:Ly
        if k>i+1
            tmp=[tmp;y{k-i}(:)-y{k-i-1}(:)];
        else
            tmp=[tmp;zeroy];
        end
    end

    u{k}=u{k-1}+rou*Fai1{k}(:,:)'*(yd{k+1}(:)-y{k}(:)-Fai{k}(:,Mu+1:My*Ly+Mu*Lu)*tmp)/(lamda+norm(Fai1{k},2)^2);
    u1(k)=u{k}(1);
    u2(k)=u{k}(2);
    %model
    y1(k+1)=2.5*y1(k)*y1(k-1)/(1+y1(k)^2+y2(k-1)^2+y1(k-2)^2)+0.7*sin(0.5*(y1(k)+y1(k-1)))*cos(0.5*(y1(k)+y1(k-1)))+.09*u1(k)*u2(k-1)+1.2*u1(k)+1.6*u1(k-2)+0.5*u2(k);
    y2(k+1)=5*y2(k)*y2(k-1)/(1+y2(k)^2+y1(k-1)^2+y2(k-2)^2)+u2(k)+1.1*u2(k-1)+1.4*u2(k-2)+0.5*u1(k);
    y{k+1}(1)=y1(k+1);
    y{k+1}(2)=y2(k+1);
end

mark=8;
step=20;
figure(1)
plot(0,'-k','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(yd1,'k-','LineWidth',2);hold on;
plot(0:k,y1,'-.b','LineWidth',2);hold on;
plot(1:step:N,y1(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;
grid on;xlabel('时刻');ylabel('跟踪性能');legend({'y_1^{*}(k)','y_1(k)'},'Interpreter','tex');
xlim([0 N]);
figure(2)
plot(0,'-k','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(yd2,'k-','LineWidth',2);hold on;
plot(0:k,y2,'-.b','LineWidth',2);hold on;
plot(1:step:N,y2(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;
grid on;xlabel('时刻');ylabel('跟踪性能');legend({'y_2^{*}(k)','y_2(k)'},'Interpreter','tex');
xlim([0 N]);
figure(3)
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'--r^','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(u1(:),'-.b','LineWidth',2);hold on;
plot(1:step:N,u1(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(u2(:),'--r','LineWidth',2);hold on;
plot(10:step:N,u2(10:step:N),'r^','MarkerSize',mark,'LineWidth',2);hold on;
grid on;xlabel('时刻');ylabel('控制输入');legend({'u_1(k)','u_2(k)'},'Interpreter','tex');
xlim([0 N]);