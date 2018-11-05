clear all;
close all;
%采样点
N=700;

%控制器参数
ny=1;
nu=2;
% % % % %projection
eita =0.2;
miu =1;
rou=0.7;
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
fai2(1,:) =[-2 0.5 0.2];
fai2(2,:)=fai2(1,:);fai2(3,:)=fai2(1,:);fai2(4,:)=fai2(1,:);fai2(5,:)=fai2(1,:);
%程序循环
for k=6:N
    if k == 8
        aaa = 1;
    end
    if k<=490
        ref(k)=(.4*(-1).^round(k/50))*0.999+0.001*ref(k-1);
        refs(k)=.4*(-1).^round(k/50);
    else
        ref(k)=(0.1+.1*(-1).^round(k/50))*0.999+0.001*ref(k-1);
        refs(k)=(0.1+0.1*(-1).^round(k/50));
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
    %fai2(6,:)=fai2(5,:);
% % % %          for the LS method
%     if k>6
%         fai2(k,:)=fai2(k-1,:)+[dy2(k-1,1:ny) du2(k-1,1:nu)]*p*(y2(k)-y2(k-1)-fai2(k-1,:)*[dy2(k-1,1:ny) du2(k-1,1:nu)]')/(alpha(k-1)+[dy2(k-1,1:ny) du2(k-1,1:nu)]*p*[dy2(k-1,1:ny) du2(k-1,1:nu)]');
%     
%         p=(p-(p*[dy2(k-1,1:ny) du2(k-1,1:nu)]'*[dy2(k-1,1:ny) du2(k-1,1:nu)]*p)/(alpha(k-1)+[dy2(k-1,1:ny) du2(k-1,1:nu)]*p*[dy2(k-1,1:ny) du2(k-1,1:nu)]'))/alpha(k-1);
%         alpha(k)=alpha0.*alpha(k-1)+(1-alpha0);
%     else
%         fai2(k,:)=fai2(k-1,:);
%         p=p;
%         alpha(k)=alpha(k-1);
%     end
    for i=1:ny
        dy2(k,i)=y2(k-i+1)-y2(k-i);
    end
    if ny<=0
        u2(k) = u2(k-1)+rou*fai2(k,ny+1)*(refs(k)-y2(k)-fai2(k,ny+2:ny+nu)*du2(k-1,2:nu)')/(lamda+fai2(k,ny+1).^2); 
    else
        u2(k) = u2(k-1)+rou*fai2(k,ny+1)*(refs(k)-y2(k)-fai2(k,1:ny)*dy2(k,:)'-fai2(k,ny+2:ny+nu)*du2(k-3,1:nu-1)')/(lamda+fai2(k,ny+1).^2); 
    end
    for i=1:nu
        du2(k,i)=u2(k-i+1)-u2(k-i);
    end
    %model
    a(k)=1;b(k)=0;
    b(k)=4*round(k./100)+sin(k/100);%
    y2(k+1)=(-0.9*a(k)*y2(k)+(b(k)+1)*u2(k))/(1+y2(k)^2);
    du2(k)=u2(k)-u2(k-1);
end

%控制器参数
nu=1;
eita =0.2;
miu =1;
rou1=0.7;
rou2=0.7;
rou3=0.7;
rou4=0.7;
rou5=0.7;
rou6=0.7;
lamda =7;
%初值
%y(1:nu+1)=0;u(1:nu)=0;du(1:nu,1:nu)=0;
ref(1:5)=0;
y(1:6)=0;y(4)=1;y(5)=0.2;
u(1:6)=0;u(5)=0.5;
% for i=1:ny
%     dy(5,i)=y(5-i+1)-y(5-i);
% end
% for i=1:nu
%     du(5,i)=u(5-i+1)-u(5-i);
% end
% y(1)=-1;y(2)=1;y(3)=0.5;y(4)=0.5;
% u(1:3)=0;
du(1:3,1:nu)=0;
% du(N,1)=0;
I=eye(nu);
%控制器伪偏导数初值
% fai(1:nu,1) =2;
% fai(1:nu,2:nu)=0;
fai(1:3,1) =[-2 0.5 0.2];
fai(1:3,2)=[-2 0.5 0.2];
fai(1:3,3)=[-2 0.5 0.2];
fai(1:3,4)=[-2 0.5 0.2];
fai(1:3,5)=[-2 0.5 0.2];
fai(1:3,6)=[-2 0.5 0.2];
xi(1:N) = 0.1;
%程序循环
%for k=nu+1:N
for k=4:N
    if k<=490
        ref(k)=(.4*(-1).^round(k/50))*0.999+0.001*ref(k-1);
        refs(k)=.4*(-1).^round(k/50);
    else
        ref(k)=(0.1+.1*(-1).^round(k/50))*0.999+0.001*ref(k-1);
        refs(k)=(0.1+0.1*(-1).^round(k/50));
    end
    a(k)=1+round(k/500);
    %H(:, k) = [y(k)-y(k-1) y(k-1)-y(k-2) du(k, 1) xi(k) xi(k - 1) xi(k - 2)]';
    H(:, k-1) = [y(k-1)-y(k-2) y(k-2)-y(k-3) du(k-1, 1) xi(k-1) xi(k - 2) xi(k - 3)]';
    xi(k) = y(k)-y(k-1) - fai(k-1, :)*H(:, k-1);
    fai(k, :)=fai(k-1, :)+(eita*(y(k)-y(k-1)- fai(k-1,:)*H(:, k-1)) .* H(:, k-1))';
    if (fai(k,1)<10^(-5)) || ((du(k-1,1:nu)*du(k-1,1:nu)')^0.5<10^(-5))
        fai(k,1)=2;
    end
    u(k) = u(k-1)+(rou3*fai(k,3)*(refs(k)-y(k)) - rou1*fai(k,1)*fai(k,3)*(y(k)-y(k-1)) - rou2*fai(k,2)*fai(k,3)*(y(k-1)-y(k-2)) - rou4*fai(k,3)*fai(k,4)*xi(k) - rou5*fai(k,3)*fai(k,5)*xi(k-1) - rou6*fai(k,3)*fai(k,6)*xi(k-2))/(lamda+fai(k,3).^2);

    for i=1:nu
        du(k,i)=u(k-i+1)-u(k-i);
    end
    %model
    a(k)=1;b(k)=0;
    b(k)=4*round(k./100)+sin(k/100);%
    y(k+1)=(-0.9*a(k)*y(k)+(b(k)+1)*u(k))/(1+y(k)^2);
    du(k)=u(k)-u(k-1);
    
%     y(k+1) = y(k+1) + rand_data(k);
end

% 求方差
var_y = sum((y(1:N) - refs).^ 2) / N
var_y2 = sum((y2(1:N) - refs).^ 2) / N

mark=8;
step=20;
figure(1)
plot(0,'-k','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(refs,'k','LineWidth',2);hold on;
plot(y,'-.b','LineWidth',2);hold on;grid on;
plot(y2,'--r','LineWidth',2);
plot(1:step:N,y(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;
grid on;xlabel('时刻');ylabel('跟踪性能');legend({'y^{*}(k)','y(k)改进','y(k)原始'},'Interpreter','tex');
xlim([0 700]);
figure(2)
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(u,'-.b','LineWidth',2);hold on;grid on;
plot(u2,'--r','LineWidth',2);hold on;grid on;
plot(1:step:N,u(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;
grid on;xlabel('时刻');ylabel('控制输入');legend({'u(k)改进','u(k)原始'});

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
ylim([-3 2.5]);
grid on;xlabel('时刻');ylabel('PG的估计值');legend('\phi_1(k)的估计值','\phi_2(k)的估计值','\phi_3(k)的估计值');