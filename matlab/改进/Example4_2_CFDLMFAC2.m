clear all;
%close all;
N=1000;

%控制器参数
nu=1;
eita =1;
miu =1;
rou1=0.6;
rou2=0.6;
lamda =2;

%初值
%y(1:nu+1)=0;u(1:nu)=0;du(1:nu,1:nu)=0;
y(1)=-1;y(2)=1;y(3)=0.5;
u(1:2)=0;
du(1:2,1:nu)=0;
%期望值
for k=1:N+1
    if k<=300
        yd(k)=5*sin(k*pi/50)+2*cos(k*pi/100);
    else
        if k<=700
            yd(k)=5*(-1)^round(k/100);
        else
            yd(k)=5*sin(k*pi/50)+2*cos(k*pi/100);
        end
     end
end
I=eye(nu);
%控制器伪偏导数初值
% fai(1:nu,1) =2;
% fai(1:nu,2:nu)=0;
fai(1:2,1) =2;
fai(1:2,2)=0;
xi(1:2) = 0;
%程序循环
for k=3:N
    a(k)=1+round(k/500);
    xi(k) = y(k)-y(k-1) - fai(k-1, :)*[du(k-1,1) xi(k - 1)]';
    fai(k, :)=fai(k-1, :)+eita*(y(k)-y(k-1)-fai(k-1,1)*du(k-1,1)' - fai(k - 1, 2) * xi(k - 1))*[du(k-1,1) xi(k - 1)];
    if (fai(k,1)<10^(-5)) || ((du(k-1,1:nu)*du(k-1,1:nu)')^0.5<10^(-5))
        fai(k,1)=2;
    end
    u(k) = u(k-1)+(rou1*fai(k,1)*(yd(k+1)-y(k)) - rou2*fai(k,1)*fai(k,2)*xi(k))/(lamda+fai(k,1).^2);
    %model
    if k<=500
        if k<=2
            y(k+1)=5*y(k)*y(k-1)/(1+y(k)^2+y(k-1)^2)+u(k)+1.2*u(k-1);
        else
            y(k+1)=5*y(k)*y(k-1)/(1+y(k)^2+y(k-1)^2+y(k-2)^2)+u(k)+1.2*u(k-1); 
        end
    else
        y(k+1)=2.5*y(k)*y(k-1)/(1+y(k)^2+y(k-1)^2)+1.2*u(k)+1.4*u(k-1)+0.7*sin(0.5*(y(k)+y(k-1)))*cos(0.5*(y(k)+y(k-1)));    
    end
    for i=1:nu
        du(k,i)=u(k-i+1)-u(k-i);
    end
    emax(k+1)=yd(k+1)-y(k+1);
end
%%%%%%%%%%%%%%%%%%%%%%%%
%控制器参数
nu=1;
eita =1;
miu =1;
rou=0.6;
lamda =2;
%初值
%y(1:nu+1)=0;u(1:nu)=0;du(1:nu,1:nu)=0;
y2(1)=-1;y2(2)=1;y2(3)=0.5;
u2(1:2)=0;
du2(1:2,1:nu)=0;
%期望值
fai2(1:2,1) =2;
fai2(1:2,2:nu)=0;
for k=3:N
    a(k)=1+round(k/500);
    fai2(k,1:nu)=fai2(k-1,1:nu)+eita*(y2(k)-y2(k-1)-fai2(k-1,1:nu)*du2(k-1,1:nu)')*du2(k-1,1:nu)/(miu+du2(k-1,1:nu)*du2(k-1,1:nu)');
    if (fai2(k,1)<10^(-5)) %|| ((du(k-1,1:nu)*du(k-1,1:nu)')^0.5<10^(-5))
        fai2(k,1)=0.5;
    end
    if nu==1
        u2(k) = u2(k-1)+rou*fai2(k,1)*(yd(k+1)-y2(k))/(lamda+fai2(k,1).^2);        
    else
        u2(k) = u2(k-1)+rou*fai2(k,1)*(yd(k+1)-y2(k)-fai2(k,2:nu)*du2(k-1,1:nu-1)')/(lamda+fai2(k,1).^2); 
    end
    %model
    if k<=500
        if k<=2
            y2(k+1)=5*y2(k)*y2(k-1)/(1+y2(k)^2+y2(k-1)^2)+u2(k)+1.2*u2(k-1);
        else
            y2(k+1)=5*y2(k)*y2(k-1)/(1+y2(k)^2+y2(k-1)^2+y2(k-2)^2)+u2(k)+1.2*u2(k-1); 
        end
    else
        y2(k+1)=2.5*y2(k)*y2(k-1)/(1+y2(k)^2+y2(k-1)^2)+1.2*u2(k)+1.4*u2(k-1)+0.7*sin(0.5*(y2(k)+y2(k-1)))*cos(0.5*(y2(k)+y2(k-1)));    
    end
    for i=1:nu
        du2(k,i)=u2(k-i+1)-u2(k-i);
    end
    emax(k+1)=yd(k+1)-y2(k+1);
end

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
grid on;xlabel('时刻');ylabel('跟踪性能');legend({'y^{*}(k)','\lambda=0.1时的输出y(k)','\lambda=2时的输出y(k)'},'Interpreter','tex');
xlim([0,1000]);ylim([-30,20]);

figure(2)
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'--r^','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(u,'-.b','LineWidth',2);hold on;
plot(1:step:N,u(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(u2,'r--','LineWidth',2);
plot(10:step:N,u2(10:step:N),'r^','MarkerSize',mark,'LineWidth',2);hold on;
grid on;xlabel('时刻');ylabel('控制输入');legend({'\lambda=0.1时的控制输入u(k)','\lambda=2时的控制输入u(k)'},'Interpreter','tex');


figure(3)
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'--r^','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(fai,'-.b','LineWidth',2);hold on;
plot(1:step:N,fai(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(fai2,'--r','LineWidth',2);grid on;
plot(10:step:N,fai2(10:step:N),'r^','MarkerSize',mark,'LineWidth',2);hold on;

xlabel('时刻');ylabel('PPD估计值');legend({'\lambda=0.1时PPD的估计值','\lambda=2时PPD的估计值'},'Interpreter','tex');
%plot(yd,'k');hold on;
%plot(0:k,y,'b');hold on;
% figure
% plot(u,'b.');hold on;
% figure
% plot(fai);hold on;

