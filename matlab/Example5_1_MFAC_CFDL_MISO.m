clear all;
%close all;
%采样时间
N=1000;

%控制器参数
nu=1;
eita =1;
miu =1;
rou=1;
lamda =3;

%初值
y(1)=1,y(2)=0.5;
u1(1)=1;du1(1:2,1:nu)=0;
u2(1)=1;du2(1:2,1:nu)=0;

epsilon=10^(-5);
M=50;

%期望值
for k=1:N+1
    %yd(k)=1;(-1)^round(k/100);
    yd(k)=(-1)^round(k/100);
end

%控制器伪偏导数初值
fai1(1:2,1) =0.5;
fai1(1:2,2:nu)=0;
fai2(1:2,1)=-0.2;
fai2(1:2,2:nu)=0;
%程序循环
for k=2:N
    fai1(k,1:nu)=fai1(k-1,1:nu)+eita*(y(k)-y(k-1)-fai1(k-1,1:nu)*du1(k-1,1:nu)'-fai2(k-1,1:nu)*du2(k-1,1:nu)')*du1(k-1,1:nu)/(miu+du1(k-1,1:nu)*du1(k-1,1:nu)'+du2(k-1,1:nu)*du2(k-1,1:nu)');
    if (abs(fai1(k,1))<epsilon)|(abs(fai1(k,1))>M)
        fai1(k,1)=fai1(1,1);
    end
    fai2(k,1:nu)=fai2(k-1,1:nu)+eita*(y(k)-y(k-1)-fai1(k-1,1:nu)*du1(k-1,1:nu)'-fai2(k-1,1:nu)*du2(k-1,1:nu)')*du2(k-1,1:nu)/(miu+du1(k-1,1:nu)*du1(k-1,1:nu)'+du2(k-1,1:nu)*du2(k-1,1:nu)');
    if (abs(fai2(k,1))<epsilon)|(abs(fai2(k,1))>M)
        fai1(k,1)=fai1(1,1);
    end
    u1(k) = u1(k-1)+rou*fai1(k,1)*(yd(k+1)-y(k)-fai1(k,2:nu)*du1(k-1,1:nu-1)'-fai2(k,2:nu)*du2(k-1,1:nu-1)')/(lamda+fai1(k,1).^2+fai2(k,1).^2); 
    u2(k) = u2(k-1)+rou*fai2(k,1)*(yd(k+1)-y(k)-fai1(k,2:nu)*du1(k-1,1:nu-1)'-fai2(k,2:nu)*du2(k-1,1:nu-1)')/(lamda+fai1(k,1).^2+fai2(k,1).^2); 

    %model
    y(k+1)=(5*y(k)+2*u1(k)-3*u2(k).^2+2*u1(k).^2)/(5+u1(k)+5*u2(k));
    
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
fai=[fai1 fai2];

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
xlim([0 1000]);
figure(2)
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'--r^','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(u1,'-.b','LineWidth',2);hold on;
plot(1:step:N,u1(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(u2,'r--','LineWidth',2);hold on;
plot(10:step:N,u2(10:step:N),'r^','MarkerSize',mark,'LineWidth',2);hold on;
grid on;xlabel('时刻');ylabel('控制输入');legend({'u_1(k)','u_2(k)'},'Interpreter','tex');
ylim([-0.5 2.2]);
figure(3)
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'--r^','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(fai(:,1),'-.b','LineWidth',2);hold on;
plot(1:step:N,fai(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(fai(:,2),'--r','LineWidth',2);grid on;
plot(10:step:N,fai(10:step:N,2),'r^','MarkerSize',mark,'LineWidth',2);hold on;
grid on;xlabel('时刻');ylabel('伪梯度估计值');legend('\phi_1(k)的估计值','\phi_2(k)的估计值');
ylim([-0.5 1.5]);
