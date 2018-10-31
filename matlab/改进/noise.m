clear all; close all;
% clc
L=500;  %仿真长度
c = [1 -0.5];
nc = length(c) - 1;
xik=zeros(nc,1);  %白噪声初值
xi=randn(L,1);  %产生均值为0，方差为1的高斯白噪声序列

for k=1:L
    e(k)=c*[xi(k);xik];  %产生有色噪声
    %数据更新
    for i=nc:-1:2
        xik(i)=xik(i-1);
    end
    xik(1)=xi(k);
end

subplot(2,1,1);
plot(xi);
xlabel('k');ylabel('噪声幅值');title('白噪声序列');
subplot(2,1,2);
plot(e);
xlabel('k');ylabel('噪声幅值');title('有色噪声序列');