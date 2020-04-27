t = linspace(0,1,1000);
fc = 10;
sig_ori = exp(1i*fc*t)*exp(1i*3);
% rng(0);
sample_time = 100;
sig = zeros(sample_time,length(sig_ori));
for j = 1:100
    for i = 1:sample_time
        tau = 1e-1*randn(1,1);
        sig(i,:) = sig_ori*exp(1i*2*pi*fc*tau);
        ang(i) = angle(sig(i,1))-angle(sig_ori(1));
    end
    ang_mean(j) = mean(ang);
end
sig_fix = sig_ori*exp(1i*-1*mean(ang));
figure;
plot(real(sig_fix));
hold on;
plot(real(sig_ori));
plot(real(sig(20,:)))
%%
sig_fix = sig(1,:)*exp(1i*-(ang(1)));
figure;plot(angle(sig(3,:)))
figure;plot(angle(sig_ori))
%%
tau = 1e-10*randn(1,100000);
test = exp(1i*2*pi*fc*tau);
ang_test = angle(test);
figure;histogram(ang_test);