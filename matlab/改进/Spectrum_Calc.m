function [Yf,f] = Spectrum_Calc(yt,Fs)
%功能：实现快速fourier变换
%输入参数：yt为时域信号序列，Fs为采样频率
%返回值：Yf为经过fft的频域序列，f为相应频率

L = length(yt);
NFFT = 2^nextpow2(L);
Yf = fft(yt,NFFT)/L;

Yf = 2*abs(Yf(1:NFFT/2+1));
f = Fs/2 * linspace(0,1,NFFT/2+1);
end