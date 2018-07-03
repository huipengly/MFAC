clear all;
%close all;
%采样时间
N=800;

%控制器参数
nu=3;
eita =0.5;
miu =1;
rou=0.5;
lamda =0.01;

%初值
y1(1:3)=0,y1(2)=1,y2(1:3)=0,y2(2)=1;
u(1:2,1:2)=0;du1(1:nu)=0;du2(1:nu)=0;du=[du1 du2];

b1=0.2;
b2=0.5;
M=b2*10;

%期望值
for k=1:N+1
    yd1(k)=5*sin(k/50)+2*cos(k/20);
    yd2(k)=2*sin(k/50)+5*cos(k/20);
end

%控制器伪偏导数初值
fai11(1:2,1) =0.5;
fai11(1:2,2:nu)=0;
fai12(1:2,1)=0;
fai12(1:2,2:nu)=0;
fai21(1:2,1)=0;
fai21(1:2,2:nu)=0;
fai22(1:2,1)=0.5;
fai22(1:2,2:nu)=0;
fai1=[fai11 fai12];
fai2=[fai21 fai22];
%程序循环
for k=3:N
%   estimate algorithm
    fai1(k,:)=fai1(k-1,:)+eita*(y1(k)-y1(k-1)-fai1(k-1,:)*du')*du/(miu+du*du');
    if (abs(fai1(k,1))<b2)|(abs(fai1(k,1))>M)|(fai1(k,1)<0)
        fai1(k,1)=fai1(1,1);
    end
    if abs(fai1(k,1+nu))>b1
        fai1(k,1+nu)=fai1(1,1+nu);
    end
    fai11(k,:)=fai1(k,1:nu);fai12(k,:)=fai1(k,nu+1:2*nu);
    
    fai2(k,:)=fai2(k-1,:)+eita*(y2(k)-y2(k-1)-fai2(k-1,:)*du')*du/(miu+du*du');
    if abs(fai2(k,1))>b1
        fai2(k,1)=fai2(1,1);
    end
    if (abs(fai2(k,1+nu))<b2)|(abs(fai2(k,1+nu))>M)|(fai2(k,1+nu)<0)
        fai2(k,1+nu)=fai2(1,1+nu);
    end
    fai21(k,:)=fai2(k,1:nu);fai22(k,:)=fai2(k,nu+1:2*nu);
%   control law
    fai=[fai11(k,1) fai12(k,1);fai21(k,1) fai22(k,1)];
    er1=yd1(k+1)-y1(k)-fai11(k,2:nu)*du1(1:nu-1)'-fai12(k,2:nu)*du2(1:nu-1)';
    er2=yd2(k+1)-y2(k)-fai21(k,2:nu)*du1(1:nu-1)'-fai22(k,2:nu)*du2(1:nu-1)';
    er=[er1 er2];
    u(k,:)=u(k-1,:)+rou*er*fai/(lamda+norm(fai,2)^2);
    %model
        y1(k+1)=2.5*y1(k)*y1(k-1)/(1+y1(k)^2+y2(k-1)^2+y1(k-2)^2)+0.7*sin(0.5*(y1(k)+y1(k-1)))*cos(0.5*(y1(k)+y1(k-1)))+.09*u(k,1)*u(k-1,2)+1.2*u(k,1)+1.6*u(k-2,1)+0.5*u(k,2);
        y2(k+1)=5*y2(k)*y2(k-1)/(1+y2(k)^2+y1(k-1)^2+y2(k-2)^2)+u(k,2)+1.1*u(k-1,2)+1.4*u(k-2,2)+0.5*u(k,1);
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

mark=8;
step=20;
figure(1)
% plot(0,'-k','MarkerSize',mark,'LineWidth',2);hold on;
% plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
% plot(0,'--r^','MarkerSize',mark,'LineWidth',2);hold on;
% set(gca,'LineWidth',2,'fontsize',28);
% plot(yd1,'k-','LineWidth',2);hold on;
% plot(0:k,y1,'-.b','LineWidth',2);hold on;
% plot(1:step:N,y1(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(0:k,y1,'--r','LineWidth',2);hold on;
plot(10:step:N,y1(10:step:N),'r^','MarkerSize',mark,'LineWidth',2);hold on;
grid on;xlabel('时刻');ylabel('跟踪性能');legend({'y_1^{*}(k)','CFDL-MFAC','PFDL-MFAC'},'Interpreter','tex');
xlim([0 N]);
figure(2)
% plot(0,'-k','MarkerSize',mark,'LineWidth',2);hold on;
% plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
% plot(0,'--r^','MarkerSize',mark,'LineWidth',2);hold on;
% set(gca,'LineWidth',2,'fontsize',28);
% plot(yd2,'k-','LineWidth',2);hold on;
% plot(0:k,y2,'-.b','LineWidth',2);hold on;
% plot(1:step:N,y2(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(0:k,y2,'--r','LineWidth',2);hold on;
plot(10:step:N,y2(10:step:N),'r^','MarkerSize',mark,'LineWidth',2);hold on;
grid on;xlabel('时刻');ylabel('跟踪性能');legend({'y_2^{*}(k)','CFDL-MFAC','PFDL-MFAC'},'Interpreter','tex');
xlim([0 N]);
figure(3)
% plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
% plot(0,'--r^','MarkerSize',mark,'LineWidth',2);hold on;
% set(gca,'LineWidth',2,'fontsize',28);
% plot(u(:,1),'-.b','LineWidth',2);hold on;
% plot(1:step:N,u(1:step:N,1),'bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(u(:,1),'--r','LineWidth',2);hold on;
plot(10:step:N,u(10:step:N,1),'r^','MarkerSize',mark,'LineWidth',2);hold on;
grid on;xlabel('时刻');ylabel('控制输入');legend({'CFDL-MFAC','PFDL-MFAC'},'Interpreter','tex');

figure(4)
% plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
% plot(0,'--r^','MarkerSize',mark,'LineWidth',2);hold on;
% set(gca,'LineWidth',2,'fontsize',28);
% plot(u(:,2),'-.b','LineWidth',2);hold on;
% plot(1:step:N,u(1:step:N,2),'bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(u(:,2),'--r','LineWidth',2);hold on;
plot(10:step:N,u(10:step:N,2),'r^','MarkerSize',mark,'LineWidth',2);hold on;
grid on;xlabel('时刻');ylabel('控制输入');legend({'CFDL-MFAC','PFDL-MFAC'},'Interpreter','tex');

