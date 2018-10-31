clear all; close all;
L=500;  %仿真长度
d=[1 -1.5 0.7 0.1]; c=[1 0.5 0.2];  % 分子分母多项式系数
nd=length(d)-1 ;nc=length(c)-1;   %阶次
xik=zeros(nc,1);  %白噪声初值
ek=zeros(nd,1);
xi=randn(L,1);  %产生均值为0，方差为1的高斯白噪声序列

for k=1:L
    e(k)=-d(2:nd+1)*ek+c*[xi(k);xik];  %产生有色噪声
    %数据更新
    for i=nd:-1:2
        ek(i)=ek(i-1);
    end
    ek(1)=e(k);
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

%测试功率谱

[y1,f1] = Spectrum_Calc(xi',512);
p1 = 1/L * y1.*conj(y1);

figure(2)
subplot(211)
plot(f1,p1)

[y2,f2] = Spectrum_Calc(e,512);
p2 = 1/L * y2.*conj(y2);
subplot(212)
plot(f2,p2)