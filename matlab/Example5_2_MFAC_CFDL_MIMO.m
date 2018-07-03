clear all;
%close all;
%采样时间
N=1500;

%控制器参数
nu=1;
eita =1;
miu =1;
rou=1;
lamda =0.5;

%初值
x11(1:2)=0.5;
x12(1:2)=0;
x21(1:2)=0.5;
x22(1:2)=0;
y1(1:2)=x11(1:2),y2(1:2)=x21(1:2);
u(1,1:2)=0;du1(1:nu)=0;du2(1:nu)=0;du=[du1 du2];

epsilon=10^(-5);
M=50;

%期望值
for k=1:N+1
    yd1(k)=0.5+0.25*cos(0.25*pi*k/100)+0.25*sin(0.5*pi*k/100);
    yd2(k)=0.5+0.25*sin(0.25*pi*k/100)+0.25*sin(0.5*pi*k/100);
end

%控制器伪偏导数初值
fai11(1:2,1) =0.5;
fai11(1:2,2:nu)=0;
fai12(1:2,1)=0.1;
fai12(1:2,2:nu)=0;
fai21(1:2,1)=0.1;
fai21(1:2,2:nu)=0;
fai22(1:2,1)=0.5;
fai22(1:2,2:nu)=0;
fai1=[fai11 fai12];
fai2=[fai21 fai22];
%程序循环
for k=2:N
    a(k)=1+0.1*sin(2*pi*k/N);
    b(k)=1+0.1*cos(2*pi*k/N);
%   estimate algorithm
    fai1(k,:)=fai1(k-1,:)+eita*(y1(k)-y1(k-1)-fai1(k-1,:)*du')*du/(miu+du*du');
    if (abs(fai1(k,1))<epsilon)|(abs(fai1(k,1))>M)
        fai1(k,1)=fai1(1,1);
    end
%     if (abs(fai1(k,1+nu))<epsilon)|(abs(fai1(k,1+nu))>M)
%         fai1(k,1+nu)=fai1(1,1+nu);
%     end
    fai11(k,:)=fai1(k,1:nu);fai12(k,:)=fai1(k,nu+1:2*nu);
    fai2(k,:)=fai2(k-1,:)+eita*(y2(k)-y2(k-1)-fai2(k-1,:)*du')*du/(miu+du*du');
    if (abs(fai2(k,2))<epsilon)|(abs(fai2(k,2))>M)
        fai2(k,1)=fai2(1,1);
    end
%     if (abs(fai2(k,1+nu))<epsilon)|(abs(fai2(k,1+nu))>M)
%         fai2(k,1+nu)=fai2(1,1+nu);
%     end
    fai21(k,:)=fai2(k,1:nu);fai22(k,:)=fai2(k,nu+1:2*nu);
%   control law
    fai=[fai11(k,1) fai12(k,1);fai21(k,1) fai22(k,1)];
    er1=yd1(k+1)-y1(k)-fai11(k,2:nu)*du1(1:nu-1)'-fai12(k,2:nu)*du2(1:nu-1)';
    er2=yd2(k+1)-y2(k)-fai21(k,2:nu)*du1(1:nu-1)'-fai22(k,2:nu)*du2(1:nu-1)';
    er=[er1 er2];
    u(k,:)=u(k-1,:)+rou*er*fai/(lamda+norm(fai,2)^2);
    %model
    x11(k+1)=x11(k)^2/(1+x11(k)^2)+0.3*x12(k);
    x12(k+1)=x11(k)^2/(1+x12(k)^2+x21(k)^2+x22(k)^2)+a(k)*u(k,1);
    x21(k+1)=x21(k)^2/(1+x21(k)^2)+0.2*x22(k);
    x22(k+1)=x21(k)^2/(1+x11(k)^2+x12(k)^2+x22(k)^2)+b(k)*u(k,2);
    y1(k+1)=x11(k+1);
    y2(k+1)=x21(k+1);
    for i=1:nu
        if (k-i)<=0
            du1(i)=0;
            du2(i)=0;
        else
            du1(i)=u(k-i+1,1)-u(k-i,1);
            du2(i)=u(k-i+1,2)-u(k-i,2);
        end
    end
    du=[du1 du2];
end
fai=[fai1 fai2];

mark=8;
step=20;
figure(1)
plot(0,'-k','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(yd1,'k-','LineWidth',2);hold on;
plot(0:k,y1,'-.b','LineWidth',2);hold on;
plot(1:step:N,y1(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;
grid on;xlabel('时刻');ylabel('y_1(k)的跟踪性能');legend({'y_1^{*}(k)','y_1(k)'},'Interpreter','tex');
xlim([0 1500]);
figure(2)
plot(0,'-k','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(yd2,'k-','LineWidth',2);hold on;
plot(0:k,y2,'-.b','LineWidth',2);hold on;
plot(1:step:N,y2(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;
grid on;xlabel('时刻');ylabel('y_2(k)的跟踪性能');legend({'y_2^{*}(k)','y_2(k)'},'Interpreter','tex');
xlim([0 1500]);
figure(3)
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'--r^','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(u(:,1),'-.b','LineWidth',2);hold on;
plot(1:step:N,u(1:step:N,1),'bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(u(:,2),'r--','LineWidth',2);hold on;
plot(10:step:N,u(10:step:N,2),'r^','MarkerSize',mark,'LineWidth',2);hold on;
grid on;xlabel('时刻');ylabel('控制输入');legend({'u_1(k)','u_2(k)'},'Interpreter','tex');

figure(4)
plot(0,'-ko','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'--r^','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,':g>','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(fai11(:,1),'-k','LineWidth',2);hold on;
plot(1:step:N,fai11(1:step:N),'ko','MarkerSize',mark,'LineWidth',2);hold on;
plot(fai12(:,1),'--r','LineWidth',2);grid on;
plot(10:step:N,fai12(10:step:N),'r^','MarkerSize',mark,'LineWidth',2);hold on;
plot(fai21(:,1),'-.b','LineWidth',2);hold on;
plot(1:step:N,fai21(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(fai22(:,1),':g','LineWidth',2);grid on;
plot(10:step:N,fai22(10:step:N),'g>','MarkerSize',mark,'LineWidth',2);hold on;
ylim([-0.1 0.8]);
grid on;xlabel('时刻');ylabel('伪雅可比矩阵估计值');legend('\phi_{11}(k)的估计值','\phi_{12}(k)的估计值','\phi_{21}(k)的估计值','\phi_{22}(k)的估计值');
 
