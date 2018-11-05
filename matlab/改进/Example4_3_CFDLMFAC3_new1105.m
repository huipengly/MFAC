clear all;
%close all;
N=1000;

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
for k=1:N+1
    yd(k)=0.5+0.5*(-1)^round(k/200);
end
for k=4:N
    %H(:, k) = [y(k)-y(k-1) y(k-1)-y(k-2) du(k, 1) xi(k) xi(k - 1) xi(k - 2)]';
    H(:, k-1) = [y(k-1)-y(k-2) y(k-2)-y(k-3) du(k-1, 1) xi(k-1) xi(k - 2) xi(k - 3)]';
    xi(k) = y(k)-y(k-1) - fai(k-1, :)*H(:, k-1);
    fai(k, :)=fai(k-1, :)+(eita*(y(k)-y(k-1)- fai(k-1,:)*H(:, k-1)) .* H(:, k-1))';
    if (fai(k,1)<10^(-5)) || ((du(k-1,1:nu)*du(k-1,1:nu)')^0.5<10^(-5))
        fai(k,1)=2;
    end
    u(k) = u(k-1)+(rou3*fai(k,3)*(yd(k)-y(k)) - rou1*fai(k,1)*fai(k,3)*(y(k)-y(k-1)) - rou2*fai(k,2)*fai(k,3)*(y(k-1)-y(k-2)) - rou4*fai(k,3)*fai(k,4)*xi(k) - rou5*fai(k,3)*fai(k,5)*xi(k-1) - rou6*fai(k,3)*fai(k,6)*xi(k-2))/(lamda+fai(k,3).^2);

    for i=1:nu
        du(k,i)=u(k-i+1)-u(k-i);
    end
    %model
    b(k)=0.1+0.1*round(k/100);
    if k<=200
        y(k+1)=1.5*y(k)-0.7*y(k-1)+0.1*(u(k)+b(k)*u(k-1));    
    else if k<=400
            y(k+1)=1.5*y(k)-0.7*y(k-1)+0.1*(u(k-2)+b(k)*u(k-3));             
        else if k<=600
                y(k+1)=1.5*y(k)-0.7*y(k-1)+0.1*(u(k-4)+b(k)*u(k-5));             
            else if k<=800
                    y(k+1)=1.5*y(k)-0.7*y(k-1)+0.1*(u(k-6)+b(k)*u(k-7));             
                else
                        y(k+1)=1.5*y(k)-0.7*y(k-1)+0.1*(u(k-8)+b(k)*u(k-9));             
                end
            end
        end
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
lamda =15;
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
    b(k)=0.1+0.1*round(k/100);
    if k<=200
        y2(k+1)=1.5*y2(k)-0.7*y2(k-1)+0.1*(u2(k)+b(k)*u2(k-1));    
    else if k<=400
            y2(k+1)=1.5*y2(k)-0.7*y2(k-1)+0.1*(u2(k-2)+b(k)*u2(k-3));             
        else if k<=600
                y2(k+1)=1.5*y2(k)-0.7*y2(k-1)+0.1*(u2(k-4)+b(k)*u2(k-5));             
            else if k<=800
                    y2(k+1)=1.5*y2(k)-0.7*y2(k-1)+0.1*(u2(k-6)+b(k)*u2(k-7));             
                else
                        y2(k+1)=1.5*y2(k)-0.7*y2(k-1)+0.1*(u2(k-8)+b(k)*u2(k-9));             
                end
            end
        end
    end
    for i=1:nu
        du2(k,i)=u2(k-i+1)-u2(k-i);
    end
    emax(k+1)=yd(k+1)-y2(k+1);
end

% 求方差
var_y = sum((y - yd).^ 2) / N
var_y2 = sum((y2 - yd).^ 2) / N

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
grid on;xlabel('时刻');ylabel('跟踪性能');legend({'y^{*}(k)','改进的输出y(k)','\lambda=15时的输出y(k)'},'Interpreter','tex');
xlim([0,1000]);ylim([-1,2.8]);

figure(2)
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'--r^','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(u,'-.b','LineWidth',2);hold on;
plot(1:step:N,u(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(u2,'r--','LineWidth',2);
plot(10:step:N,u2(10:step:N),'r^','MarkerSize',mark,'LineWidth',2);hold on;
ylim([-1,2.7]);
grid on;xlabel('时刻');ylabel('控制输入');legend({'改进的控制输入u(k)','\lambda=15时的控制输入u(k)'},'Interpreter','tex');


figure(3)
plot(0,'-.bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(0,'--r^','MarkerSize',mark,'LineWidth',2);hold on;
set(gca,'LineWidth',2,'fontsize',28);
plot(fai,'-.b','LineWidth',2);hold on;
plot(1:step:N,fai(1:step:N),'bs','MarkerSize',mark,'LineWidth',2);hold on;
plot(fai2,'--r','LineWidth',2);grid on;
plot(10:step:N,fai2(10:step:N),'r^','MarkerSize',mark,'LineWidth',2);hold on;
ylim([-1,2.5]);
xlabel('时刻');ylabel('PPD估计值');legend({'改进的PPD的估计值','\lambda=15时PPD的估计值'},'Interpreter','tex');


%plot(yd,'k');hold on;
%plot(0:k,y,'b');hold on;
% figure
% plot(u,'b.');hold on;
% figure
% plot(fai);hold on;

