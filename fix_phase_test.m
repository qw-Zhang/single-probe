fc = 2.535e9;
t = linspace(0,1e-9,100);
% doppler frequency fd = fc*v*cos(theta)/c
v = 10; %m/s
c = 3e8;
theta = 30*pi/180;
fd = fc*v*cos(theta)/c;
sample_time = 1000;
sig = zeros(sample_time,length(t));
p = linspace(-pi,pi,100);
ph = (-pi + 2*pi*rand(1));
ph
sig_ori = exp(1i*2*pi*(fc + fd)*t)*exp(1i*ph);
sig_ori_fft = fft(sig_ori,512);
[v_ori_max, pos_ori_max] = max(abs(sig_ori_fft));
ang_ori = angle(sig_ori_fft(pos_ori_max));
tau = 1e-11*randn(1,sample_time);
ang = angle(exp(1i*2*pi*fc*tau));
[ang_est,delta_ang] = deal(zeros(1,sample_time));
for i = 1:sample_time
    sig(i,:) = sig_ori*exp(1i*2*pi*fc*tau(i));
    fft_temp = fft(sig(i,:),512);
    [v_max, pos_max] = max(abs(fft_temp));
    ang_est(i) = angle(fft_temp(pos_max));
%     delta_ang(i) = ang_est(i) - ang_ori;
%     if delta_ang(i) > pi
%         delta_ang(i) = 2*pi - delta_ang(i);
%     end
%     if delta_ang(i) < -pi
%         delta_ang(i) = 2*pi + delta_ang(i);
%     end
end
figure;histogram(ang_est);
tau_est = delta_ang/(2*pi*fc);
ang_mean = mean(ang);
ang_est_mean = mean(ang_est);
[h_y,h_x] = hist(ang_est,50);
[vh_y,ph_y] = max(h_y);
[v_min,pos_min] = min(abs(ang_est - mean(ang_est)));
sig_res = sig(1,:)*exp(1i*-ang_est(1))*exp(1i*h_x(ph_y));
% sig_res = sig(pos_min,:)*exp(1i*-v_min);
figure;
plot(imag(sig_res)); hold on; plot(imag(sig_ori));