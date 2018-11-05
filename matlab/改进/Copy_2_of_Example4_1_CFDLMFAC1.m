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
    e(k)=-d(2:nd+1)*ek+c*[xii(k);xik];  %产生有色噪声
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
rand_data = zeros(N, 1);
% 使用随机噪声
% rand_data = rand(N, 1);
% 使用白噪声
% rand_data = xii;
% 使用有色噪声
% rand_data = e;

% 噪声减小10倍
rand_data = rand_data ./ 100;

%控制器参数
nu=1;
eita =0.2;
miu =1;
rou1=0.6;
rou2=0.6;
rou3=0.6;
rou4=0.6;
rou5=0.6;
rou6=0.6;
lamda =7;

%初值
%y(1:nu+1)=0;u(1:nu)=0;du(1:nu,1:nu)=0;
y(1)=-1;y(2)=1;y(3)=0.5;y(4)=0.5;
u(1:3)=0;
du(1:3,1:nu)=0;
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
fai(1:3,1) =2;
fai(1:3,2)=0.1;
fai(1:3,3)=0.1;
fai(1:3,4)=0.1;
fai(1:3,5)=0.1;
fai(1:3,6)=0.1;
xi(1:3) = 0.1;
%程序循环
%for k=nu+1:N
for k=4:N
    a(k)=1+round(k/500);
    H(:, k) = [y(k)-y(k-1) y(k-1)-y(k-2) du(k-1, 1) xi(k - 1) xi(k - 2) xi(k - 3)]';
    xi(k) = y(k)-y(k-1) - fai(k-1, :)*H(:, k);
    fai(k, :)=fai(k-1, :)+(eita*(y(k)-y(k-1)-fai(k-1,1)*du(k-1,1)' - fai(k-1,:)*H(:, k)) .* H(:, k))';
    if (fai(k,1)<10^(-5)) || ((du(k-1,1:nu)*du(k-1,1:nu)')^0.5<10^(-5))
        fai(k,1)=2;
    end
    u(k) = u(k-1)+(rou3*fai(k,3)*(yd(k+1)-y(k)) - rou1*fai(k,1)*fai(k,3)*(y(k)-y(k-1)) - rou2*fai(k,2)*fai(k,3)*(y(k-1)-y(k-2)) - rou4*fai(k,3)*fai(k,4)*xi(k) - rou5*fai(k,3)*fai(k,5)*xi(k-1) - rou6*fai(k,3)*fai(k,6)*xi(k-2))/(lamda+fai(k,3).^2);
    %model
    if k<=500
        y(k+1)=y(k)/(1+y(k).^2)+u(k)^3;
    else
        y(k+1)=(y(k)*y(k-1)*y(k-2)*u(k-1)*(y(k-2)-1)+a(k)*u(k))/(1+y(k-1)^2+y(k-2)^2);    
    end
    for i=1:nu
        du(k,i)=u(k-i+1)-u(k-i);
    end
    y(k+1) = y(k+1) + rand_data(k);
    emax(k+1)=yd(k+1)-y(k+1);
end

%控制器参数
ny=2;
nu=1;
% % % % %projection
eita =0.2;
miu =1;
rou=0.6;
lamda =7;
% % % % % %LS
% p=1000*eye(3);alpha(5)=.95;alpha0=.99;
%初值
ref(1:5)=0;
y2(1:6)=0;y2(4)=1;y2(5)=0.2;
u2(1:6)=0;u2(5)=0.5;
for i=1:ny
    dy2(5,i)=y2(5-i+1)-y2(5-i);
end
for i=1:nu
    du2(5,i)=u2(5-i+1)-u2(5-i);
end
I=eye(nu);
%控制器伪偏导数初值
fai2(1,:) =[2 0.5 0.2];
fai2(2,:)=fai2(1,:);fai2(3,:)=fai2(1,:);fai2(4,:)=fai2(1,:);fai2(5,:)=fai2(1,:);
%程序循环
%for k=nu+1:N
for k=3:N
    if k<=490
        ref(k)=(5*(-1).^round(k/100))*0.45+0.55*y(k);
        refs(k)=5*(-1).^round(k/100);
    else
        refs(k)=3.5+0.5*(-1).^round(k/100);
        ref(k)=refs(k)*0.45+0.55*y(k);
    end
% % % %          for the projection method
    if ny<=0
        fai2(k,:)=fai2(k-1,:)+eita*(y2(k)-y2(k-1)-fai2(k-1,:)*[du2(k-1,1:nu)]')*[du2(k-1,1:nu)]/(miu+[du2(k-1,1:nu)]*[du2(k-1,1:nu)]');
    else
        fai2(k,:)=fai2(k-1,:)+eita*(y2(k)-y2(k-1)-fai2(k-1,:)*[dy2(k-1,1:ny) du2(k-1,1:nu)]')*[dy2(k-1,1:ny) du2(k-1,1:nu)]/(miu+[dy2(k-1,1:ny) du2(k-1,1:nu)]*[dy2(k-1,1:ny) du2(k-1,1:nu)]');
    end
    if fai2(k,1+ny)<0.00001
        fai2(k,1+ny)=0.5;
    end
    fai2(6,:)=fai2(5,:);
    for i=1:ny
        dy2(k,i)=y2(k-i+1)-y2(k-i);
    end
    if ny<=0
        u2(k) = u2(k-1)+rou*fai2(k,ny+1)*(refs(k)-y2(k)-fai2(k,ny+2:ny+nu)*du2(k-1,2:nu)')/(lamda+fai2(k,ny+1).^2); 
    else
        u2(k) = u2(k-1)+rou*fai2(k,ny+1)*(refs(k)-y2(k)-fai2(k,1:ny)*dy2(k,:)'-fai2(k,ny+2:ny+nu)*du2(k-1,1:nu-1)')/(lamda+fai2(k,ny+1).^2); 
    end
    for i=1:nu
        du2(k,i)=u2(k-i+1)-u2(k-i);
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
grid on;xlabel('时刻');ylabel('跟踪性能');legend({'y^{*}(k)','改进方法','\lambda=0.1时的输出y(k)'},'Interpreter','tex');
xlim([0,1000]);

figure(2)
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'--r^','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(u,'-.b','LineWidth',2);hold on;
plot(1:step:N,u(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(u2,'r--','LineWidth',2);
plot(10:step:N,u2(10:step:N),'r^','MarkerSize',mark,'LineWidth',2);hold on;
grid on;xlabel('时刻');ylabel('控制输入');legend({'改进方法','\lambda=0.1时的控制输入u(k)'},'Interpreter','tex');


figure(3)
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'--r^','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(fai,'-.b','LineWidth',2);hold on;
plot(1:step:N,fai(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(fai2,'--r','LineWidth',2);grid on;
plot(10:step:N,fai2(10:step:N),'r^','MarkerSize',mark,'LineWidth',2);hold on;
xlabel('时刻');ylabel('PPD估计值');legend({'改进方法','\lambda=0.1时PPD的估计值'},'Interpreter','tex');

 