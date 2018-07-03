clear all;
%close all;
%采样时间
h=0.001;N=400;T=1;

%控制器参数
Kp=0.15;Ti=0.5;Td=0;

%初值
y(1:6)=0;y(4)=1;
u(1:5)=0;
e(1:6)=5;
ee=0;

%期望值
for k=1:N+1
    yd(k)=5*(-1)^round(k/80);
end

%程序循环
for k=6:N
    ee=ee+e(k);
    u(k) = Kp*(e(k)+T*ee/Ti);
    %model
    if k<=200
        y(k+1)=2.5*y(k)*y(k-1)/(1+y(k).^2+y(k-1).^2)+0.7*sin(0.5*(y(k)+y(k-1)))+1.4*u(k-1)+1.2*u(k);
    else
        y(k+1)=-0.1*y(k)-0.2*y(k-1)-0.3*y(k-2)+0.1*u(k)+0.02*u(k-1)+0.03*u(k-2);    
    end
    e(k+1)=yd(k+1)-y(k+1);
end
mark=8;
step=20;
figure(1)
plot(0,'-k','MarkerSize',mark','LineWidth',2);hold on;
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'--r^','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(0:k,yd,'k','LineWidth',2);hold on;
plot(0:k,y,'-.b','LineWidth',2);hold on;
plot(1:step:N,y(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;

figure(2)
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'--r^','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(u,'-.b','LineWidth',2);hold on;
plot(1:step:N,u(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;
