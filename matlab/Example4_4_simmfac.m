clear all;
%close all;
%采样时间
h=0.001;N=400;

%控制器参数
nu=3;
eita =0.75;
miu =0.9;
rou=0.5;
lamda =0.01;

%初值
y(1:6)=0;y(4)=1;
u(1:5)=0;
du(1:5,1:nu)=0;

%期望值
for k=1:N+1
    yd(k)=5*(-1)^round(k/80);
end
%控制器伪偏导数初值
fai(1:5,1) =1;
fai(1:5,2:nu)=0;
%程序循环
for k=6:N
    fai(k,1:nu)=fai(k-1,1:nu)+eita*(y(k)-y(k-1)-fai(k-1,1:nu)*du(k-1,1:nu)')*du(k-1,1:nu)/(miu+du(k-1,1:nu)*du(k-1,1:nu)');
    if fai(k,1)<0.00001
        fai(k,1)=0.5;
    end
    if nu==1
         u(k) = u(k-1)+rou*fai(k,1)*(yd(k+1)-y(k))/(lamda+fai(k,1).^2);       
    else
        u(k) = u(k-1)+rou*fai(k,1)*(yd(k+1)-y(k)-fai(k,2:nu)*du(k-1,1:nu-1)')/(lamda+fai(k,1).^2); 
    end
    %model
    if k<=200
        y(k+1)=2.5*y(k)*y(k-1)/(1+y(k).^2+y(k-1).^2)+0.7*sin(0.5*(y(k)+y(k-1)))*cos(0.5*(y(k)+y(k-1)))+1.4*u(k-1)+1.2*u(k);
    else
        y(k+1)=-0.1*y(k)-0.2*y(k-1)-0.3*y(k-2)+0.1*u(k)+0.02*u(k-1)+0.03*u(k-2);    
    end
    for i=1:nu
        du(k,i)=u(k-i+1)-u(k-i);
    end
    emax(k+1)=yd(k+1)-y(k+1);
end
mark=8;
step=20;
% plot(yd,'k');hold on;
figure(1)
plot(0:k,y,'--r','LineWidth',2);hold on;
plot(10:step:N,y(10:step:N),'r^','MarkerSize',mark,'LineWidth',2);hold on;
ylim([-8 10]);
grid on;xlabel('时刻');ylabel('跟踪性能');legend({'y^{*}(k)','PID','PFDL-MFAC'},'Interpreter','tex');
figure(2)
plot(u,'r-.','LineWidth',2);hold on;
plot(10:step:N,u(10:step:N),'r^','MarkerSize',mark,'LineWidth',2);hold on;
ylim([-60 80]);
grid on;xlabel('时刻');ylabel('控制输入');legend({'PID','PFDL-MFAC'},'Interpreter','tex');
figure(3)
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'--r^','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,':g>','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(0:k-1,fai(:,1),'-.b','LineWidth',2);
plot(5:step:N,fai(5:step:N,1),'bs','MarkerSize',mark,'LineWidth',2);hold on;

plot(0:k-1,fai(:,2),'--r','LineWidth',2);
plot(10:step:N,fai(10:step:N,2),'r^','MarkerSize',mark,'LineWidth',2);hold on;

plot(0:k-1,fai(:,3),':g','LineWidth',2);
plot(15:step:N,fai(15:step:N,3),'g>','MarkerSize',mark,'LineWidth',2);hold on;

grid on;xlabel('时刻');ylabel('PG的估计值');legend('\phi_1(k)的估计值','\phi_2(k)的估计值','\phi_3(k)的估计值');

