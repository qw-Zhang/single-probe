% plot_2ray_exp_model.m

%MIMO-OFDM Wireless Communications with MATLAB  Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
%?2010 John Wiley & Sons (Asia) Pte Ltd

clear;
scale=1e-9;                         % ns
Ts=10*scale;                        % Sampling time
t_rms=10*scale;                     % RMS delay spread
num_ch=10000;                       % # of channel
%%
% 2-ray model
pow_2=[0.5 0.5];  delay_2=[0 t_rms*2]/scale;
H_2 = Ray_model(num_ch).'*sqrt(pow_2);
avg_pow_h_2 = mean(H_2.*conj(H_2));
subplot(121), stem(delay_2,pow_2), hold on, stem(delay_2,avg_pow_h_2,'r.');
xlabel('Delay[ns]'), ylabel('Channel Power[linear]');
title('Ideal PDP and simulated PDP of 2-ray model');
legend('Ideal','Simulation');  axis([0 140 0 0.7]);
%%
% Exponential model
pow_e=exp_pdp(t_rms,Ts);
% delay_e=(0:length(pow_e)-1)*Ts/scale;
delay_e = [0,5,7,12,15,20];
pow_db = [-3,-4.3,-5.7,-7.3,-9,-11.4];  
pow_e = [1, 0.7, 0.6, 0.4, 0.3, 0.15];
H_e = Ray_model(num_ch).'*sqrt(pow_e);
avg_pow_h_e = mean(H_e.*conj(H_e));
H = zeros(num_ch,length(pow_e));
for i = 1:6
    ray_h(i,:) = Ray_model(num_ch);
    temp_H = ray_h(i,:)*sqrt(pow_e(i));
    H((delay_e(i)+1):num_ch,i) = temp_H(1:(num_ch - delay_e(i)));
    % i need test the signal convolved H and generate the synchornation
    % peaks. And signal is filtered by awgn channel.
end
Hf = sum(H,2); %still that question...one column data to caculate PDP...
H_awgn = awgn(Hf,0);
figure;hold on;
for i = 1:6
    [a,b] = xcorr(H_awgn,ray_h(i,:));    
    plot(b,abs(a));
    v(i) = max(abs(a));
end
v = v./v(1);
v = v.^2
pow_H = mean(H.*conj(H));
figure;
stem(avg_pow_h_e);hold on;stem(pow_H);stem(v);
% avg_pow_h_e = mean(power(abs(H_e),2));
% subplot(122), stem(delay_e,pow_e), hold on, stem(delay_e,avg_pow_h_e,'r.');
% xlabel('Delay[ns]'), ylabel('Channel Power[linear]');
% title('Ideal PDP and simulated PDP of exponential model');
% legend('Ideal','Simulation');  
% axis([0 140 0 0.7]);

P = 4;				% modulation order
hMod   = comm.QPSKModulator;
hDemod = comm.QPSKDemodulator('OutputDataType','double');
s = RandStream.create('mt19937ar', 'seed',55408);
prevStream = RandStream.setGlobalStream(s);
frmLen = 2048;
data = randi([0 P-1], frmLen, 1);

% Modulate data
modData = step(hMod, data); 
% corr_test = zeros(length(delay_e),length(modData));
for i = 1:length(delay_e)
    tx(i,:) = delayseq(conv(modData, H_e(:,1)*sqrt(pow_e(i))), delay_e(i));
    corr_test(i,:) = sqrt(pow_e(i))*delayseq(awgn(modData,0), delay_e(i));
end
% tx_sig = sum(tx,1);
tx_sig = sum(corr_test,1);

avg_pow_sig = mean(tx'.*conj(tx'));
avg_pow_sig = avg_pow_sig*(avg_pow_h_e(1)/avg_pow_sig(1));

%%
%gcc广义互相关
% % N=1024;  %����
% % Fs=500;  %����Ƶ��
% % n=0:N-1;
% % t=n/Fs;   %ʱ������
% % d=2;     %�ӳٵ���
% % x1=5*cos(2*pi*10*n/Fs);     %�ź�1
% % x2=5*cos(2*pi*10*(n+2)/Fs); %�ź�2
% % X1=fft(x1,2*N-1);
X1 = fft(modData,frmLen*2-1);
% X2=fft(x2,2*N-1);
X2 = fft(tx_sig,frmLen*2-1);
G12=X1.*conj(X2);
G11=X1.*conj(X1);
G22=X2.*conj(X2);
Pxy=G12./sqrt(G11.*G22);
Cxy=fftshift(ifft(Pxy));
% Rl2=fix(length(Cxy)/2);
% t1=(-Rl2:Rl2)/Fs;
plot(abs(Cxy(:,1)))

X1 = fft(modData,frmLen*2-1);
X2 = fft(tx_sig,frmLen*2-1);
X1J=conj(X1);
X2J=conj(X2);
X12=X2.*X1J;
R12=(ifft(X12));

                 
% Rxy=Pxy./(X1.*conj(X1)+eps);                                          %1/G11, ROTH��Ȩ
% % Rxy=Pxy./(sqrt((X1.*conj(X1)).*(X2.*conj(X2)))+eps);         %1/sqrt(G11*G2), SCOT��Ȩ
% % Rxy=Pxy./(abs(Pxy)+eps);                                                 %1/abs(G12), PHAT��Ȩ
% 
% R2=ifft(Rxy);
% R2=real(R2);
% 
% R2=fftshift(R2);               
% %ȡ���ֵ
% c2=max(R2);
% a2=find(R2==c2) ;

% Rl2=fix(length(R2)/2);
% delay2=(a2-(Rl2)-1)/fs  ;%������źŵ��ӳ�ʱ��
% 
% %PATH��Ȩ
% R1=G12./abs(G12);
% %SCOT��Ȩ��Ȩ
% Pxy=G12./sqrt(G11.*G22);
%%
a = xcorr(tx_sig,modData);
% figure;plot(b,real(a));
h_temp = a(length(tx_sig)-1:length(tx_sig)+999);
% for i = i:length(delay_e)
%     h(i) =  H_e(:,1)*sqrt(pow_e(i))./h_temp';
% end

[h_a, h_b] = xcorr(h_temp, H_e(:,1));
figure;stem(h_b,abs(h_a));
%%
for i = 6:1
    temp(i) = power(cursor_info(i).Position(2),2);
end
temp/temp(6)
%  power(temp,2)
%%
N = 2;M = 2;
H = (randn(frmLen, N, M) + 1i*randn(frmLen, N, M))/sqrt(2);
for i = 1:2
    for j = 1:2
        h_pdp(:,:,i,j) = H(:,i,j)*sqrt(pow_e);
        avg_pow_h(:,i,j) = mean(h_pdp(:,:,i,j).*conj(h_pdp(:,:,i,j)));
    end
end

for i = 1:2
    for j = 1:2
        avg_pow_h_22(:,i,j) = avg_pow_h(:,i,j)/avg_pow_h(1,i,j);
    end
end
