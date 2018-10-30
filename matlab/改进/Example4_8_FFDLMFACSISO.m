clear all;
close all;
%采样点
N=500;

% 新添加参数
% p = 0.1;
% t0 = -p^3;
% t1 = 3*p^2;
% t2 = -3*p;
% t3 = 1;
t0 = -0.25;
t1 = 0.75;
t2 = -1;
t3 = 1;

a0 = 1;

%控制器参数
ny=2;
nu=1;
% % % % %projection
eita =0.2;
miu =1;
rou=0.7;
% rho_init = [4825.95464138029;16289.8329435659;-9478.60068589456];
rho_init = [0.7 0.7 0.7];
rho = [1.29441544976025;-36.8007707938611;-2.58883089952050]
rho1 = rho_init(1);
rho2 = rho_init(2);
rho3 = rho_init(3);

lamda =1;
% % % % % %LS
% p=1000*eye(3);alpha(5)=.95;alpha0=.99;
%初值
% ref(1:5)=0;
y(1:6)=0;y(4)=1;y(5)=0.2;
u(1:6)=0;u(5)=0.5;
for i=1:ny
    dy(5,i)=y(5-i+1)-y(5-i);
end
for i=1:nu
    du(5,i)=u(5-i+1)-u(5-i);
end
I=eye(nu);
%控制器伪偏导数初值
fai(1,:) =[2 0.5 0.2];
fai(2,:)=fai(1,:);fai(3,:)=fai(1,:);fai(4,:)=fai(1,:);fai(5,:)=fai(1,:);
%程序循环
for k=6:N
    if (k == 300)
        rho1 = rho(1);
        rho2 = rho(2);
        rho3 = rho(3);
    end
    
    refs(k) = 0.5;
%     if k<=490
% %         ref(k)=(5*(-1).^round(k/100))*0.45+0.55*y(k);
%         refs(k)=5*(-1).^round(k/100);
%     else
%         refs(k)=3.5+0.5*(-1).^round(k/100);
% %         ref(k)=refs(k)*0.45+0.55*y(k);
%     end
% % % %          for the projection method
    if ny<=0
        fai(k,:)=fai(k-1,:)+eita*(y(k)-y(k-1)-fai(k-1,:)*[du(k-1,1:nu)]')*[du(k-1,1:nu)]/(miu+[du(k-1,1:nu)]*[du(k-1,1:nu)]');
    else
        fai(k,:)=fai(k-1,:)+eita*(y(k)-y(k-1)-fai(k-1,:)*[dy(k-1,1:ny) du(k-1,1:nu)]')*[dy(k-1,1:ny) du(k-1,1:nu)]/(miu+[dy(k-1,1:ny) du(k-1,1:nu)]*[dy(k-1,1:ny) du(k-1,1:nu)]');
    end
    if fai(k,1+ny)<0.00001
        fai(k,1+ny)=0.5;
    end
    fai(6,:)=fai(5,:);
    for i=1:ny
        dy(k,i)=y(k-i+1)-y(k-i);
    end
%     if ny<=0
%         u(k) = u(k-1)+rou*fai(k,ny+1)*(refs(k)-y(k)-fai(k,ny+2:ny+nu)*du(k-1,2:nu)')/(lamda+fai(k,ny+1).^2); 
%     else
        u(k) = u(k-1)+(rho3*fai(k,ny+1)*(refs(k)-y(k))-fai(k,ny+1)*[rho1 rho2].*fai(k,1:ny)*dy(k,:)')/(lamda+fai(k,ny+1).^2); 
%     end
    for i=1:nu
        du(k,i)=u(k-i+1)-u(k-i);
    end
    
    %model
    a(k)=1;b(k)=0;
    b(k)=4*round(k./100)+sin(k/100);%
    y(k+1)=(-0.9*a(k)*y(k)+(b(k)+1)*u(k))/(1+y(k)^2);
    du(k)=u(k)-u(k-1);
    
%     
%     % 新添加的
%     phi1 = fai(k, 1);
%     phi2 = fai(k, 2);
%     phi3 = fai(k, 3);
% 
%     a1 = -phi1;
%     a2 = -phi2;
%     b1 =  phi3;
%     f0(k) = t0 / a0;
%     g0 = (t1 - f0(k) * (a1 - a0)) / b1;
%     g1 = (t2 - f0(k) * (a2 - a1)) / b1;
%     g2 = (t3 + f0(k) * a2) / b1;
%     h1k = zeros(3, 3);
%     temp1 = (lamda + phi3 ^ 2);
%     h1k(1, 1) = f0(k) * phi1 * phi3 / temp1;
%     h1k(1, 3) = f0(k) * phi3 / temp1;
%     h1k(2, 1) = -f0(k) * phi1 * phi3 / temp1;
%     h1k(2, 2) = f0(k) * phi2 * phi3 / temp1;
%     h1k(3, 2) = -f0(k) * phi2 * phi3 / temp1;
%     h1(:, :, k) = h1k;
%     rho(:, :, k) = h1k \ [g0 g1 g2]';        % inv(h1k) * [g0 g1 g2]'
% 
%     rho1 = rho(1, k);
%     rho2 = rho(2, k);
%     rho3 = rho(3, k);

end


% 新添加的
phi1 = fai(k, 1);
phi2 = fai(k, 2);
phi3 = fai(k, 3);

a1 = -phi1;
a2 = -phi2;
b1 =  phi3;
f0(k) = t0 / a0;
g0 = (t1 - f0(k) * (a1 - a0)) / b1;
g1 = (t2 - f0(k) * (a2 - a1)) / b1;
g2 = (t3 + f0(k) * a2) / b1;
h1k = zeros(3, 3);
temp1 = (lamda + phi3 ^ 2);
h1k(1, 1) = f0(k) * phi1 * phi3 / temp1;
h1k(1, 3) = f0(k) * phi3 / temp1;
h1k(2, 1) = -f0(k) * phi1 * phi3 / temp1;
h1k(2, 2) = f0(k) * phi2 * phi3 / temp1;
h1k(3, 2) = -f0(k) * phi2 * phi3 / temp1;
h1(:, :, k) = h1k;
rho = h1k \ [g0 g1 g2]';        % inv(h1k) * [g0 g1 g2]'

% figure
mark=8;
step=20;
figure(1)
plot(0,'-k','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(refs,'k','LineWidth',2);hold on;
plot(y,'-.b','LineWidth',2);hold on;grid on;
plot(1:step:N,y(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;
grid on;xlabel('时刻');ylabel('跟踪性能');legend({'y^{*}(k)','y(k)'},'Interpreter','tex');
xlim([0 700]);
figure(2)
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(u,'-.b','LineWidth',2);hold on;grid on;
plot(1:step:N,u(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;
grid on;xlabel('时刻');ylabel('控制输入');

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
ylim([-0.5,2.2]);
grid on;xlabel('时刻');ylabel('PG的估计值');legend('\phi_1(k)的估计值','\phi_2(k)的估计值','\phi_3(k)的估计值');