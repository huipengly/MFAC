clear all;
%close all;
%采样点
N=700;

%控制器参数
ny=1;
nu=2;
% % % % %projection
eita =0.2;
miu =1;
rou=0.7;
lamda =0.1;
% % % % % %LS
% p=1000*eye(3);alpha(5)=.95;alpha0=.99;
%初值
ref(1:5)=0;
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
    if k<=490
        ref(k)=(5*(-1).^round(k/100))*0.45+0.55*y(k);
        refs(k)=5*(-1).^round(k/100);
    else
        refs(k)=3.5+0.5*(-1).^round(k/100);
        ref(k)=refs(k)*0.45+0.55*y(k);
    end
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
% % % %          for the LS method
%     if k>6
%         fai(k,:)=fai(k-1,:)+[dy(k-1,1:ny) du(k-1,1:nu)]*p*(y(k)-y(k-1)-fai(k-1,:)*[dy(k-1,1:ny) du(k-1,1:nu)]')/(alpha(k-1)+[dy(k-1,1:ny) du(k-1,1:nu)]*p*[dy(k-1,1:ny) du(k-1,1:nu)]');
%     
%         p=(p-(p*[dy(k-1,1:ny) du(k-1,1:nu)]'*[dy(k-1,1:ny) du(k-1,1:nu)]*p)/(alpha(k-1)+[dy(k-1,1:ny) du(k-1,1:nu)]*p*[dy(k-1,1:ny) du(k-1,1:nu)]'))/alpha(k-1);
%         alpha(k)=alpha0.*alpha(k-1)+(1-alpha0);
%     else
%         fai(k,:)=fai(k-1,:);
%         p=p;
%         alpha(k)=alpha(k-1);
%     end
    for i=1:ny
        dy(k,i)=y(k-i+1)-y(k-i);
    end
    if ny<=0
        u(k) = u(k-1)+rou*fai(k,ny+1)*(refs(k)-y(k)-fai(k,ny+2:ny+nu)*du(k-1,2:nu)')/(lamda+fai(k,ny+1).^2); 
    else
        u(k) = u(k-1)+rou*fai(k,ny+1)*(refs(k)-y(k)-fai(k,1:ny)*dy(k,:)'-fai(k,ny+2:ny+nu)*du(k-1,1:nu-1)')/(lamda+fai(k,ny+1).^2); 
    end
    for i=1:nu
        du(k,i)=u(k-i+1)-u(k-i);
    end
    %model
    if abs(u(k))<=0.2
        y(k+1)=1.6*y(k)-0.63*y(k-1);
    else
        y(k+1)=1.6*y(k)-.63*y(k-1)+(u(k)-.2*sign(u(k)))-0.5*(u(k-1)-.2*sign(u(k-1)));
    end
    du(k)=u(k)-u(k-1);
end
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