clear all; clc; close all;
N=1024;  %长度
Fs=500;  %采样频率
n=0:N-1;
t=n/Fs;   %时间序列
a1=5;     %信号幅度
a2=5;
d=9;     %延迟点数
x1=a1*cos(2*pi*10*n/Fs);     %信号1
x1=awgn(x1,20);
x2=a2*cos(2*pi*10*(n+d)/Fs); %信号2
x2=awgn(x2,20);
subplot(211);
plot(t,x1,'r');
axis([-0.2 2 -6 6]);
hold on;
plot(t,x2,':');
axis([-0.2 2 -6 6]);
legend('x1信号', 'x2信号');
xlabel('时间/s');ylabel('x1(t) x2(t)');
title('原始信号');grid on;
hold off
%互相关函数
X1=fft(x1,2*N-1);
X2=fft(x2,2*N-1);
Sxy=X1.*conj(X2);
Cxy=fftshift(ifft(Sxy));     %基本互相关
Cxy=fftshift(real(ifft(Sxy./abs(Sxy))));
subplot(212);
t1=((0:2*N-2)-N+1)/Fs;
plot(t1,Cxy,'b');
title('Rx1x2');xlabel('t/s');ylabel('Rx1x2(t)');grid on
[max0,location]=max(Cxy);     
d=location-N
Delay=d/Fs              %求得时间延迟