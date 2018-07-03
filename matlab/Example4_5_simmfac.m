clear all;
%close all;
%采样时间
N=800;

%控制器参数
nu=3;
eita =0.75;
miu =0.9;
rou=0.4;
lamda =.01;

%初值
y(1:6)=0;y(4)=1;
u(1:5)=0;
du(1:5,1:nu)=0;

%期望值
for k=1:N+1
    yd(k)=0.5*(-1)^round(k/50);
end
I=eye(nu);
%控制器伪偏导数初值
fai(1:5,1) =1;
fai(1:5,3:nu)=0;
%程序循环
for k=6:N
    fai(k,1:nu)=fai(k-1,1:nu)+eita*(y(k)-y(k-1)-fai(k-1,1:nu)*du(k-1,1:nu)')*du(k-1,1:nu)/(miu+du(k-1,1:nu)*du(k-1,1:nu)');
    if fai(k,1)<0.00001
        fai(k,1)=0.5;
    end
    u(k) = u(k-1)+rou*fai(k,1)*(yd(k+1)-y(k)-fai(k,2:nu)*du(k-1,1:nu-1)')/(lamda+fai(k,1).^2); 

    %model
    if k<=100
        a(k)=0;
    elseif k<=300
        a(k)=1;
    elseif k<=500
        a(k)=2;
    elseif k<=800
        a(k)=3;
    end
%     a(k)=1;
    y(k+1)=(y(k)*y(k-1)*y(k-2)*u(k-1)*(y(k)-1)+(1+a(k))*u(k))/(1+y(k)^2+y(k-1)^2+y(k-2)^2);
    for i=1:nu
        du(k,i)=u(k-i+1)-u(k-i);
    end
    emax(k+1)=yd(k+1)-y(k+1);
end
mark=8;
step=20;
figure(1)
plot(0,'-k','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(0:k,yd,'k','LineWidth',2);hold on;
plot(0:k,y,'-.b','LineWidth',2);hold on;
plot(1:step:N,y(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;
grid on;xlabel('时刻');ylabel('跟踪性能');legend({'y^{*}(k)','y(k)'},'Interpreter','tex');

figure(2)
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(u,'-.b','LineWidth',2);hold on;
plot(1:step:N,u(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;
grid on;xlabel('时刻');ylabel('控制输入');
% 
% plot(fai);

