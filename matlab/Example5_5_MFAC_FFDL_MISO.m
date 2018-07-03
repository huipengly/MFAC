clear all;
%close all;
%采样时间
N=1000;

%控制器参数
ny=1;
nu=1;
eita =1;
miu =1;
rou=1;
lamda =22;

%初值
y(1)=.1,y(2)=1;
u1(1)=.1;du1(1:2,1:nu)=0;
u2(1)=.2;du2(1:2,1:nu)=0;

epsilon=10^(-5);
M=50;

%期望值
for k=1:N+1
    yd(k)=0.5*(-1)^round(k/100);
    a(k)=round(k/50);
    b(k)=sin(k/100);
end

%控制器伪偏导数初值
faiy(1:2,1)=1;faiy(1:2,2:ny)=0;
fai1(1:2,1) =-.5;fai1(1:2,2)=-.5;
fai1(1:2,3:nu)=0;
fai2(1:2,1)=1;fai2(1:2,2)=0;
fai2(1:2,3:nu)=0;
%程序循环
for k=2:N
    for i=1:ny
        if (k-i)<=0
            dy(k,i)=0;
            dy(k,i)=0;
        else
            dy(k,i)=y(k-i+1)-y(k-i);
            dy(k,i)=y(k-i+1)-y(k-i);
        end
    end
    
    faiy(k,1:ny)=faiy(k-1,1:ny)+eita*(y(k)-y(k-1)-faiy(k-1,1:ny)*dy(k-1,1:ny)'-fai1(k-1,1:nu)*du1(k-1,1:nu)'-fai2(k-1,1:nu)*du2(k-1,1:nu)')*dy(k-1,1:ny)/(miu+du1(k-1,1:nu)*du1(k-1,1:nu)'+du2(k-1,1:nu)*du2(k-1,1:nu)'+dy(k-1,1:ny)*dy(k-1,1:ny)');
    
    fai1(k,1:nu)=fai1(k-1,1:nu)+eita*(y(k)-y(k-1)-faiy(k-1,1:ny)*dy(k-1,1:ny)'-fai1(k-1,1:nu)*du1(k-1,1:nu)'-fai2(k-1,1:nu)*du2(k-1,1:nu)')*du1(k-1,1:nu)/(miu+du1(k-1,1:nu)*du1(k-1,1:nu)'+du2(k-1,1:nu)*du2(k-1,1:nu)'+dy(k-1,1:ny)*dy(k-1,1:ny)');
    if (abs(fai1(k,1))<epsilon)|(abs(fai1(k,1))>M)
        fai1(k,1)=fai1(1,1)
    end
    
    fai2(k,1:nu)=fai2(k-1,1:nu)+eita*(y(k)-y(k-1)-faiy(k-1,1:ny)*dy(k-1,1:ny)'-fai1(k-1,1:nu)*du1(k-1,1:nu)'-fai2(k-1,1:nu)*du2(k-1,1:nu)')*du2(k-1,1:nu)/(miu+du1(k-1,1:nu)*du1(k-1,1:nu)'+du2(k-1,1:nu)*du2(k-1,1:nu)'+dy(k-1,1:ny)*dy(k-1,1:ny)');
    if (abs(fai2(k,1))<epsilon)|(abs(fai2(k,1))>M)
        fai2(k,1)=fai2(1,1);
    end
    
    u1(k) = u1(k-1)+rou*fai1(k,1)*(yd(k+1)-y(k)-faiy(k,1:ny)*dy(k,1:ny)'-fai1(k,2:nu)*du1(k-1,1:nu-1)'-fai2(k,2:nu)*du2(k-1,1:nu-1)')/(lamda+fai1(k,1).^2+fai2(k,1).^2); 
    u2(k) = u2(k-1)+rou*fai2(k,1)*(yd(k+1)-y(k)-faiy(k,1:ny)*dy(k,1:ny)'-fai1(k,2:nu)*du1(k-1,1:nu-1)'-fai2(k,2:nu)*du2(k-1,1:nu-1)')/(lamda+fai1(k,1).^2+fai2(k,1).^2); 

    %model
    y(k+1)=y(k)*u1(k).^2+a(k)*u2(k)-b(k)*u1(k-1)*y(k-1).^2+u2(k-1).^2;
   
    for i=1:nu
        if (k-i)<=0
            du1(k,i)=0;
            du2(k,i)=0;
        else
            du1(k,i)=u1(k-i+1)-u1(k-i);
            du2(k,i)=u2(k-i+1)-u2(k-i);
        end
    end
    emax(k+1)=yd(k+1)-y(k+1);
end
fai=[faiy fai1(:,1),fai2(:,1)];
mark=8;
step=20;
figure(1)
plot(0,'-k','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(yd,'k-','LineWidth',2);hold on;
plot(0:k,y,'-.b','LineWidth',2);hold on;
plot(1:step:N,y(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;
grid on;xlabel('时刻');ylabel('跟踪性能');legend({'y^{*}(k)','y(k)'},'Interpreter','tex');
xlim([0 N]);
figure(2)
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'--r^','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(u1,'-.b','LineWidth',2);hold on;
plot(1:step:N,u1(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(u2,'r--','LineWidth',2);hold on;
plot(10:step:N,u2(10:step:N),'r^','MarkerSize',mark,'LineWidth',2);hold on;
grid on;xlabel('时刻');ylabel('控制输入');legend({'u_1(k)','u_2(k)'},'Interpreter','tex');

figure(3)
plot(0,'-ko','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'--r^','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(fai(:,1),'-k','LineWidth',2);hold on;
plot(1:step:N,fai(1:step:N),'ko','MarkerSize',mark,'LineWidth',2);hold on;
plot(fai(:,2),'--r','LineWidth',2);grid on;
plot(10:step:N,fai(10:step:N,2),'r^','MarkerSize',mark,'LineWidth',2);hold on;
plot(fai(:,3),'-.b','LineWidth',2);hold on;
plot(1:step:N,fai(1:step:N,3),'bs','MarkerSize',mark,'LineWidth',2);hold on;
ylim([-1 2.5]);
grid on;xlabel('时刻');ylabel('伪分块梯度估计值');legend('\phi_{1}(k)的估计值','\phi_{12}(k)的估计值','\phi_{22}(k)的估计值');
 