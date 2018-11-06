% 修改phi的方法
clear all;
close all;
N=1000;

% 产生噪声 xii是白噪声，e是有色噪声
d=[1 -1.5 0.7 0.1]; c=[1 0.5 0.2];  % 分子分母多项式系数
nd=length(d)-1 ;nc=length(c)-1;   %阶次
xik=zeros(nc,1);  %白噪声初值
ek=zeros(nd,1);
xii=randn(N,1);  %产生均值为0，方差为1的高斯白噪声序列

for k=1:N
    e(k)=-d(2:nd+1)*ek+c*[xii(k);xik+0.5];  %产生有色噪声
    %数据更新
    for i=nd:-1:2
        ek(i)=ek(i-1);
    end
    ek(1)=e(k);
    for i=nc:-1:2
        xik(i)=xik(i-1);
    end
    xik(1)=xii(k);
end
% 无噪声
% rand_data = zeros(N, 1);
% 使用随机噪声
% rand_data = rand(N, 1);
% 使用白噪声
% rand_data = xii;
% 使用有色噪声
rand_data = e;

% 噪声减小100倍
rand_data = rand_data ./ 100;

%控制器参数
nu=1;
eita =1;
miu =1;
rou=0.6;
lamda =0.1;

%初值
%y(1:nu+1)=0;u(1:nu)=0;du(1:nu,1:nu)=0;
y(1)=-1;y(2)=1;y(3)=0.5;
u(1:2)=0;
du(1:2,1:nu)=0;
%期望值
for k=1:N+1
    if k<=300
        yd(k)=0.5*(-1)^round(k/100);
    else
        if k<=700
            yd(k)=0.5*sin(k*pi/100)+0.3*cos(k*pi/50);
        else
            yd(k)=0.5*(-1)^round(k/100);
        end
     end
end
I=eye(nu);
%控制器伪偏导数初值
% fai(1:nu,1) =2;
% fai(1:nu,2:nu)=0;
fai(1:2,1) =2;
fai(1:2,2:nu)=0;
%程序循环
%for k=nu+1:N
for k=3:N
    a(k)=1+round(k/500);
    fai_act(k - 1) = (y(k)-y(k-1)) / du(k-1,1:nu);
    fai(k,1:nu)=fai_act(k-1)+eita*(y(k)-y(k-1)-fai_act(k-1)*du(k-1,1:nu)')*du(k-1,1:nu)/(miu+du(k-1,1:nu)*du(k-1,1:nu)');
    if (fai(k,1)<10^(-5)) || ((du(k-1,1:nu)*du(k-1,1:nu)')^0.5<10^(-5))
        fai(k,1)=2;
    end
    if nu==1
        u(k) = u(k-1)+rou*fai(k,1)*(yd(k+1)-y(k))/(lamda+fai(k,1).^2);        
    else
        u(k) = u(k-1)+rou*fai(k,1)*(yd(k+1)-y(k)-fai(k,2:nu)*du(k-1,1:nu-1)')/(lamda+fai(k,1).^2); 
    end
    %model
    if k<=500
        y(k+1)=y(k)/(1+y(k).^2)+u(k)^3;
    else
        y(k+1)=(y(k)*y(k-1)*y(k-2)*u(k-1)*(y(k-2)-1)+a(k)*u(k))/(1+y(k-1)^2+y(k-2)^2);    
    end
    y(k+1) = y(k+1) + rand_data(k);
    for i=1:nu
        du(k,i)=u(k-i+1)-u(k-i);
    end
    emax(k+1)=yd(k+1)-y(k+1);
end

%控制器参数
nu=1;
eita =1;
miu =1;
rou=0.6;
lamda =0.1;

%初值
%y(1:nu+1)=0;u(1:nu)=0;du(1:nu,1:nu)=0;
y2(1)=-1;y2(2)=1;y2(3)=0.5;
u2(1:2)=0;
du2(1:2,1:nu)=0;

I=eye(nu);
%控制器伪偏导数初值
% fai(1:nu,1) =2;
% fai(1:nu,2:nu)=0;
fai2(1:2,1) =2;
fai2(1:2,2:nu)=0;
%程序循环
%for k=nu+1:N
for k=3:N
    a(k)=1+round(k/500);
    fai2(k,1:nu)=fai2(k-1,1:nu)+eita*(y2(k)-y2(k-1)-fai2(k-1,1:nu)*du2(k-1,1:nu)')*du2(k-1,1:nu)/(miu+du2(k-1,1:nu)*du2(k-1,1:nu)');
    if (fai2(k,1)<10^(-5)) || ((du2(k-1,1:nu)*du2(k-1,1:nu)')^0.5<10^(-5))
        fai2(k,1)=2;
    end
    if nu==1
        u2(k) = u2(k-1)+rou*fai2(k,1)*(yd(k+1)-y2(k))/(lamda+fai2(k,1).^2);        
    else
        u2(k) = u2(k-1)+rou*fai2(k,1)*(yd(k+1)-y2(k)-fai2(k,2:nu)*du2(k-1,1:nu-1)')/(lamda+fai2(k,1).^2); 
    end
    %model
    if k<=500
        y2(k+1)=y2(k)/(1+y2(k).^2)+u2(k)^3;
    else
        y2(k+1)=(y2(k)*y2(k-1)*y2(k-2)*u2(k-1)*(y2(k-2)-1)+a(k)*u2(k))/(1+y2(k-1)^2+y2(k-2)^2);    
    end
    y2(k+1) = y2(k+1) + rand_data(k);
    for i=1:nu
        du2(k,i)=u2(k-i+1)-u2(k-i);
    end
    emax2(k+1)=yd(k+1)-y2(k+1);
end

% 求方差
var_y = sum((y - yd).^ 2) / N
var_y2 = sum((y2 - yd).^ 2) / N

err_y = cumsum((y - yd).^ 2);
err_y2 = cumsum((y2 - yd).^ 2);

clf;
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
grid on;xlabel('时刻');ylabel('跟踪性能');legend({'y^{*}(k)','改进：\lambda=0.1时的输出y(k)','\lambda=0.1时的输出y(k)'},'Interpreter','tex');
xlim([0,1000]);

figure(2)
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'--r^','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(u,'-.b','LineWidth',2);hold on;
plot(1:step:N,u(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(u2,'r--','LineWidth',2);
plot(10:step:N,u2(10:step:N),'r^','MarkerSize',mark,'LineWidth',2);hold on;
grid on;xlabel('时刻');ylabel('控制输入');legend({'改进：\lambda=0.1时的控制输入u(k)','\lambda=0.1时的控制输入u(k)'},'Interpreter','tex');


figure(3)
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'--r^','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(fai,'-.b','LineWidth',2);hold on;
plot(1:step:N,fai(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(fai2,'--r','LineWidth',2);grid on;
plot(10:step:N,fai2(10:step:N),'r^','MarkerSize',mark,'LineWidth',2);hold on;
xlabel('时刻');ylabel('PPD估计值');legend({'改进：\lambda=0.1时PPD的估计值','\lambda=0.1时PPD的估计值'},'Interpreter','tex');

figure(4)
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'--r^','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(err_y,'-.b','LineWidth',2);hold on;
plot(1:step:N,err_y(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(err_y2,'--r','LineWidth',2);grid on;
plot(10:step:N,err_y2(10:step:N),'r^','MarkerSize',mark,'LineWidth',2);hold on;
xlabel('时刻');ylabel('方差和');legend({'改进方法','原始'},'Interpreter','tex');

figure()
plot(rand_data)

