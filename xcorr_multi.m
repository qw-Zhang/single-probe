%%
Txdata = load('D:\study\OTA\expriment_iecas\DATA\DATA0730\Tx_data\Txdata.dat');
temp_ch_pow = '10';
angle_map = [1 -13 146 -30 -11 -1];
%name(i,:) = ['D:\study\OTA\expriment_iecas\DATA\DATA0730\turn\data',temp_ch,temp_ch_pow,'.dat'];
%     name(i,:) = ['D:\tripod\','1',temp_ch_pow,'dbm','.dat'];
name = 'D:\study\OTA\expriment_iecas\DATA\DATA0810\multi_tripod\data1-15dbm.dat';
%     var_name = ['data',temp_ch,'10dbm'];
load(name);

%%
load('D:\study\OTA\expriment_iecas\DATA\DATA0810\multi_tripod\Txdata_3VST.dat');
for i = 1:3
    probe(i,:) = Txdata_3VST(2*i-1,:)+1i*Txdata_3VST(2*i,:);
end
%%
power_map = [-15 -10 0 5 10];
figure;
for i = 1:3
    for j = 1:3
        name_p = ['D:\study\OTA\expriment_iecas\DATA\DATA0810\multi_tripod\p',num2str((i)),'_',num2str(power_map(j)),'dbm.dat'];
        p_vna(:,(5*j-4):(5*j)) = load(name_p);
        s21(:,j) = p_vna(:,(5*j-4+1)) + 1i*p_vna(:,(5*j-4+2));
        s31(:,j) = p_vna(:,(5*j-4+3)) + 1i*p_vna(:,(5*j-4+4));
        ifft_s21(:,j) = ifft(s21(:,j));
        ifft_s31(:,j) = ifft(s31(:,j));
    end
    mean_ifft_s21(:,i) = mean(ifft_s21,2);
    mean_ifft_s31(:,i) = mean(ifft_s31,2);
    max_s21(i) = max(abs(mean_ifft_s21(:,i)));
    max_s31(i) = max(abs(mean_ifft_s31(:,i)));
    plot(abs(mean_ifft_s31(:,i)));
    hold on;
end
%%
for i = 1:3
   probe_cal_1(i,:) = probe(i,:).*max_s31(i); 
   probe_cal_2(i,:) = probe(i,:).*max_s21(i); 
end
probe_syn(1,:) = sum(probe_cal_1,1);
probe_syn(2,:) = sum(probe_cal_2,1);
%%
%     temp(1,:) = X110dbm(1,:);
%     temp(2,: )= X110dbm(2,:);
new_data = data1_15dbm(1,:) + 1i*data1_15dbm(2,:);
com_tx = Txdata(1,:)+1i*Txdata(2,:);
res_cor_com = xcorr(new_data,com_tx);
%     temp(1,:) = temp(1,:).*w';
%      r1(i,:) = xcorr(temp(1,:)*(v(1,1)/v(1,i)),Txdata(1,:));
%     r1(i,:) = xcorr(temp(1,:),Txdata(1,:));
%     [v_max,p_max] = max(abs(r1(i,:)));
%     if(r1(i,p_max)<0)
%         r1(i,:) = -1*r1(i,:);
%         new_data_re(i,:) = -1*new_data_re(i,:);
%     end
%     res = resample(r(i,:),10,1);
temp_p = int32(length(res_cor_com)/2-10);
%     plot(r1(i,temp_p:(temp_p+100)));
figure;
plot(abs(res_cor_com));

%%
res_cor_com = xcorr(probe_syn(1,:),com_tx);
temp_p = int32(length(res_cor_com)/2-10);
figure;
plot(abs(res_cor_com));

res_cor_com = xcorr(probe_syn(2,:),com_tx);
temp_p = int32(length(res_cor_com)/2-10);
figure;
plot(abs(res_cor_com));

%%
    res_new = res_cor_com(2764801-50:2764830+50);
    la = linspace(-50/30.72e6,(130-51)/30.72e6,130);
    hold on;
    
    res_new_interp = resample(abs(res_new),10000,3072);
    res_new_interp_1 = resample(abs(res_new_interp),10,1);
    
        normal = max(abs(res_new_interp));
    res_new_interp = abs(res_new_interp)./normal;
    
    la_interp = linspace((1-165)/100,(424-165)/100,424);
    la_interp_1 = linspace((1-1629)/1e9,(4240-1629)/1e9,4240);
    
    
    %     plot(la_interp, abs(res_new_interp(i,:)));
    %     plot(la_interp,  abs(res_new_interp(i ,:)))

%%
% figure;
res_new(7,:) = res_cor_com(2764801-50:2764830+50);
res_new_interp(7,:) = resample(abs(res_new(7,:)),10000,3072);
res_new_interp(7,:) = res_new_interp(7,:)./max(res_new_interp(7,:));
res_new_interp_1(7,:) = resample(abs(res_new_interp(7,:)),10,1);
plot(la_interp,  abs(res_new_interp(7 ,:)),'LineWidth',5);hold on;
xlabel('Delay [microsecond]');ylabel('Normalization Amplitude')
%%
rr_std = zeros(1,length(la_interp));
% delay_std = [0 2.05e-7 2.85e-7 6.6e-7 8.05e-7 9.25e-7];
delay_std_1 = [0 2.0e-7 2.8e-7 6.6e-7 8.0e-7 9.2e-7];
res_std = [1 0.5346 0.7471 0.3718 0.2504 0.1435];
delay = delay_std_1*1e8 + 165;
for i = 1:6
    rr_std(delay(i)) = res_std(i);
end
stem(la_interp,rr_std);
axis([-0.5 1.5 0.05 1.05 ]);
legend('Synthesized data with correct','Model standard');