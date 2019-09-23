%generate the chirp signal and qpsk signal
f0 = 10;f1 = 2000;
t = linspace(0,9999,10000);
chirp = sin(2 * pi*(f0*t / 10000 + (f1 / 2)*power(t / 10000, 2))) + 1i*sin(2 * pi*(f0*t / 10000 + (f1 / 2)*power(t / 10000, 2)));
data=[0 1 0 1 1 1 0 0 1 1]; % information
%Number_of_bit=1024;
%data=randint(Number_of_bit,1);
data_NZR=2*data-1; % Data Represented at NZR form for QPSK modulation
s_p_data=reshape(data_NZR,2,length(data)/2);  % S/P convertion of data
br=10.^6; %Let us transmission bit rate  1000000
f=br; % minimum carrier frequency
T=1/br; % bit duration
t=T/100:T/100:T; % Time vector for one bit information
%QPSK modulatio
y=[];
y_in=[];
y_qd=[];
for i=1:length(data)/2
    y1=s_p_data(1,i)*cos(2*pi*f*t); % inphase component
    y2=s_p_data(2,i)*sin(2*pi*f*t) ;% Quadrature component
    y_in=[y_in y1]; % inphase signal vector
    y_qd=[y_qd y2]; %quadrature signal vector
%     y=[y y1+y2]; % modulated signal vector
end
qpsk = y_in + 1i*y_qd;
temp = zeros(1,100);
signal = [chirp,temp,qpsk];
% signal = [chirp,qpsk];
signal = awgn(signal,100);

%%
rate = 20;
signal_in = resample(signal,rate,1);
signal_ex_len = length(signal_in)/rate;

piont = 8;

j = 1;
for i = 1:signal_ex_len
    signal_piont(j,1) = signal_in(rate*(i-1)+piont);
    j = j + 1;
end

signal_piont_in = resample(signal_piont,rate,1);
for r = 1:rate
    j = 1;
    for i = 1:signal_ex_len
        signal_piont_ex(j,r) = signal_piont_in(rate*(i-1)+r);
        j = j + 1;
    end
    [a,b] = xcorr(signal_piont_ex(:,r),signal);
    [v_in_m(2,r),p_temp] = max(abs(a));
    v_in_m(1,r) = b(p_temp);
end
% figure;plot(v_in_m(2,:));
[r_m(4),r_m(3)] = max(v_in_m(2,:));
%%
piont_move = 7;
tau = (piont_move/2000)*2*pi;
% ss = signal_piont_in*exp(1i*tau);
ss = delayseq(signal_piont_in,piont_move);
j = 1;
for i = 1:signal_ex_len
    signal_result_ex(j) = ss(rate*(i-1) + 1);
    j = j + 1;
end

% signal_result_ex = [0,signal_result_ex];
figure;
hold on;
plot(imag(signal_piont));
plot(imag(signal_result_ex));
plot(imag(signal));
axis([2000 2050 -1.5 1.5]);
title(num2str(piont_move));
legend('signal piont','signal result','signal')