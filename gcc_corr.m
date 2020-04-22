clear all; clc; close all;
N=1024;  %����
Fs=500;  %����Ƶ��
n=0:N-1;
t=n/Fs;   %ʱ������
a1=5;     %�źŷ���
a2=5;
d=9;     %�ӳٵ���
x1=a1*cos(2*pi*10*n/Fs);     %�ź�1
x1=awgn(x1,20);
x2=a2*cos(2*pi*10*(n+d)/Fs); %�ź�2
x2=awgn(x2,20);
subplot(211);
plot(t,x1,'r');
axis([-0.2 2 -6 6]);
hold on;
plot(t,x2,':');
axis([-0.2 2 -6 6]);
legend('x1�ź�', 'x2�ź�');
xlabel('ʱ��/s');ylabel('x1(t) x2(t)');
title('ԭʼ�ź�');grid on;
hold off
%����غ���
X1=fft(x1,2*N-1);
X2=fft(x2,2*N-1);
Sxy=X1.*conj(X2);
Cxy=fftshift(ifft(Sxy));     %���������
Cxy=fftshift(real(ifft(Sxy./abs(Sxy))));
subplot(212);
t1=((0:2*N-2)-N+1)/Fs;
plot(t1,Cxy,'b');
title('Rx1x2');xlabel('t/s');ylabel('Rx1x2(t)');grid on
[max0,location]=max(Cxy);     
d=location-N
Delay=d/Fs              %���ʱ���ӳ�