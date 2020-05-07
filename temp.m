fc = 2.535e9;
c = 3e8;
lambda = c/fc;
t = linspace(0,1e-9,100);
% doppler frequency fd = fc*v*cos(theta)/c
v = 10; %m/s
theta = 30*pi/180;
fd = fc*v*cos(theta)/c;
sig = exp(1i*2*pi*(fc + fd)*t);
sig_ori_fft = fft(sig,512);
[v_ori_max, pos_ori_max] = max(abs(sig_ori_fft));
ang_ori = angle(sig_ori_fft(pos_ori_max))

I = sig;
for i = 1:1000
In = awgn(sig,5);
sig_ori_fft = fft(In,512);
[v_ori_max, pos_ori_max] = max(abs(sig_ori_fft));
ang(i) = angle(sig_ori_fft(pos_ori_max));
end
histogram(ang);
% 计算信噪比函数
% I :original signal
% In:noisy signal(ie. original signal + noise signal)
snr=0;
Ps=sum(sum((I-mean(mean(I))).^2));%signal power
Pn=sum(sum((I-In).^2));           %noise power
snr=10*log10(Ps/Pn);

%%
o = open('stat_1e-10.fig');
lh = findall(o,'type','line');
y1 = lh.YData;
o1 = open('stat_1e-10_2.fig');
lh1 = findall(o1,'type','line');
y2 = lh1.YData;
figure;hold on; plot([8,16,36,72],y2);plot([8,16,36,72],y1);
xlabel('number of locations');
ylabel('RMSE');
legend('original','corrected');