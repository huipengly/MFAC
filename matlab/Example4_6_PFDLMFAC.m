clear all;
%close all;
%采样时间
h=0.001;N=800;

%控制器参数
nu=1;
eita =0.5;
miu =1;
rou=0.5;
lamda =0.01;

%初值
y(1:6)=0;y(4)=1;
u(1:5)=0;
du(1:5,1:nu)=0;

%期望值
for k=1:N+1
    yd(k)=5*sin(k/50)+2*cos(k/20);
end
I=eye(nu);
%控制器伪偏导数初值
fai(1:5,2:nu)=0;
fai(1:5,1) =1;
fai(1:5,2)=0;
fai(1:5,3)=0;
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
    if k<=400
        y(k+1)=2.5*y(k)*y(k-1)/(1+y(k)^2+y(k-1)^2)+0.7*sin(0.5*(y(k)+y(k-1)))*cos(0.5*(y(k)+y(k-1)))+.09*u(k)*u(k-1)+1.2*u(k)+1.6*u(k-2);
    else
        y(k+1)=5*y(k)*y(k-1)/(1+y(k)^2+y(k-1)^2+y(k-2)^2)+u(k)+1.1*u(k-1);
    end
    for i=1:nu
        du(k,i)=u(k-i+1)-u(k-i);
    end
    emax(k+1)=yd(k+1)-y(k+1);
end

%控制器参数
nu=3;
eita =0.5;
miu =1;
rou=0.5;
lamda =0.01;

%初值
y2(1:6)=0;y2(4)=1;
u2(1:5)=0;
du2(1:5,1:nu)=0;

%期望值
for k=1:N+1
    yd(k)=5*sin(k/50)+2*cos(k/20);
end
I=eye(nu);
%控制器伪偏导数初值
fai2(1:5,2:nu)=0;
fai2(1:5,1) =1;
fai2(1:5,2)=0;
fai2(1:5,3)=0;
%程序循环
for k=6:N
    fai2(k,1:nu)=fai2(k-1,1:nu)+eita*(y2(k)-y2(k-1)-fai2(k-1,1:nu)*du2(k-1,1:nu)')*du2(k-1,1:nu)/(miu+du2(k-1,1:nu)*du2(k-1,1:nu)');
    if fai2(k,1)<0.00001
        fai2(k,1)=0.5;
    end
    if nu==1
        u2(k) = u2(k-1)+rou*fai2(k,1)*(yd(k+1)-y2(k))/(lamda+fai2(k,1).^2);        
    else
        u2(k) = u2(k-1)+rou*fai2(k,1)*(yd(k+1)-y2(k)-fai2(k,2:nu)*du2(k-1,1:nu-1)')/(lamda+fai2(k,1).^2); 
    end
    %model
    if k<=400
        y2(k+1)=2.5*y2(k)*y2(k-1)/(1+y2(k)^2+y2(k-1)^2)+0.7*sin(0.5*(y2(k)+y2(k-1)))*cos(0.5*(y2(k)+y2(k-1)))+.09*u2(k)*u2(k-1)+1.2*u2(k)+1.6*u2(k-2);
    else
        y2(k+1)=5*y2(k)*y2(k-1)/(1+y2(k)^2+y2(k-1)^2+y2(k-2)^2)+u2(k)+1.1*u2(k-1);
    end
    for i=1:nu
        du2(k,i)=u2(k-i+1)-u2(k-i);
    end
    emax2(k+1)=yd(k+1)-y2(k+1);
end

% figure(1)
% plot(yd,'k-');hold on;
% plot(0:k,y,'b:');hold on;
% plot(0:k,y2,'r--');
% grid on;xlabel('Step');ylabel('Tracking performance');legend({'$\quad y^{*}(k)$','$\quad y(k)\quad with\quad CFDL-MFAC$','$\quad y(k)\quad with\quad PFDL-MFAC$'},'Interpreter','latex');
% figure(2)
% plot(u,'b:');hold on;plot(u2,'r--');
% grid on;xlabel('Step');ylabel('Control input');legend({'$\quad u(k)\quad with\quad CFDL-MFAC$','$\quad u(k)\quad with \quad PFDL-MFAC$'},'Interpreter','latex');
% figure(3)
% plot(fai(:,1),'b');hold on;
% plot(fai2);grid on;
% xlabel('Step');ylabel('PPD estimated value');legend({'$\quad \hat\phi_c(k)\quad in \quad CFDL-MFAC$','$\quad \hat\phi_1(k)\quad in \quad PFDL-MFAC$','$\quad \hat\phi_2(k)\quad in \quad PFDL-MFAC$','$\quad \hat\phi_3(k)\quad in \quad PFDL-MFAC$'},'Interpreter','latex');

mark=8;
step=20;
figure(1)
plot(0,'-k','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'--r^','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(yd,'k-','LineWidth',2);hold on;
plot(0:k,y,'-.b','LineWidth',2);hold on;
plot(1:step:N,y(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(0:k,y2,'--r','LineWidth',2);
plot(10:step:N,y2(10:step:N),'r^','MarkerSize',mark,'LineWidth',2);hold on;
xlim([0 800]);
grid on;xlabel('时刻');ylabel('跟踪性能');legend({'y^{*}(k)','CFDL-MFAC','PFDL-MFAC'},'Interpreter','tex');


figure(2)
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'--r^','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(u,'-.b','LineWidth',2);hold on;
plot(1:step:N,u(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(u2,'r--','LineWidth',2);
plot(10:step:N,u2(10:step:N),'r^','MarkerSize',mark,'LineWidth',2);hold on;
grid on;xlabel('时刻');ylabel('控制输入');legend({'CFDL-MFAC','PFDL-MFAC'},'Interpreter','tex');


figure(3)
plot(0,'-ko','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'--r^','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,':g>','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(0:k-1,fai(:,1),'-k','LineWidth',2);
plot(1:step:N,fai(1:step:N,1),'ko','MarkerSize',mark,'LineWidth',2);hold on;

plot(0:k-1,fai2(:,1),'-.b','LineWidth',2);
plot(5:step:N,fai2(5:step:N,1),'bs','MarkerSize',mark,'LineWidth',2);hold on;

plot(0:k-1,fai2(:,2),'--r','LineWidth',2);
plot(10:step:N,fai2(10:step:N,2),'r^','MarkerSize',mark,'LineWidth',2);hold on;

plot(0:k-1,fai2(:,3),':g','LineWidth',2);
plot(15:step:N,fai2(15:step:N,3),'g>','MarkerSize',mark,'LineWidth',2);hold on;grid on;
ylim([-2 13]);
xlabel('时刻');ylabel('PPD和PG的估计值');legend('\phi_c(k)的估计值','\phi_1(k)的估计值','\phi_2(k)的估计值','\phi_3(k)的估计值');
