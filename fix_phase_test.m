
fc = 2.535e9;
t = linspace(0,1e-9,100);
% doppler frequency fd = fc*v*cos(theta)/c
v = 10; %m/s
c = 3e8;
theta = 30*pi/180;
fd = fc*v*cos(theta)/c;

% rng(0);
sample_time = 1000;
sig = zeros(sample_time,length(t));
p = linspace(-pi,pi,100);
delta = zeros(100,100);
ph = (-pi + 2*pi*rand(1));
sig_ori = exp(1i*2*pi*(fc + fd)*t)*exp(1i*ph);
sig_ori_fft = fft(sig_ori,512);
[v_ori_max, pos_ori_max] = max(abs(sig_ori_fft));
ang_ori = angle(sig_ori_fft(pos_ori_max));
for i = 1:sample_time
    tau(i) = 1e-10*randn(1,1);
    sig(i,:) = sig_ori*exp(1i*2*pi*fc*tau(i));
    ang(i) = angle(exp(1i*2*pi*fc*tau(i)));
    fft_temp = fft(sig(i,:),512);
    [v_max, pos_max] = max(abs(fft_temp));
    ang_est(i) = angle(fft_temp(pos_max));
    delta_ang(i) = (ang_est(i) - ang_ori);
    tau_est(i) = delta_ang(i)/(2*pi*fc);
end
ang_mean = mean(ang);
ang_est_mean = mean(ang_est);
[v_min,pos_min] = min(ang_est - mean(ang_est));
sig_res = sig(pos_min,:)*exp(1i*-(ang_est_mean));
plot(real(sig_res));
hold on;
plot(real(sig_ori));