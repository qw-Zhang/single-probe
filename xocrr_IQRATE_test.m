%%
clear;
%%
%loading data
clear name var_name;
Txdata = load('D:\study\OTA\expriment_iecas\DATA\DATA0730\Tx_data\Txdata.dat');
% dir_dic = 'D:\study\OTA\expriment_iecas\DATA\DATA0821';
dir_dic = 'D:\study\OTA\expriment_iecas\DATA\DATA0822';
IQ_rate = 10;
for i = 1:6
    temp_ch = num2str(i);
    temp_IQ_rate = num2str(IQ_rate);
    name(i,:) = [dir_dic,'\data',temp_ch,'_chirp',temp_IQ_rate,'M.dat'];
    var_name(i,:) = ['data',temp_ch,'_chirp',temp_IQ_rate,'M'];
    load(name(i,:));
end

f0 = 10;f1 = 2000;
t = linspace(1,10000,10000);
chirp = sin(2 * pi*(f0*t / 10000 + (f1 / 2)*power(t / 10000, 2)));
%%
chirp_r = resample(chirp,9000,3072);
new_data_re = zeros(6,length(eval(var_name(1,:))));
% chirp_r = resample(chirp,10000,3072);
temp = eval(var_name(3,:));
new_data(3,:) =(temp(3,:) + 1i*temp(4,:));
[find,lags] = xcorr(real(new_data(3,:)),chirp_r);
figure;plot(lags,find)
len = length(real(new_data(2,:)));
[value_max,p_max] = max(abs(find(1:3/2*len)));
pp = lags(p_max)
%%
%xcorr data  of change IQrate
clear res_cor_com;
new_data = zeros(6,length(eval(var_name(1,:))));
chirp_r = resample(chirp,100*IQ_rate,3072);
com_tx = Txdata(1,:) + 1i*Txdata(2,:);
com_tx_r = resample(com_tx,100*IQ_rate,3072);
% w = hamming(length(data110dbm(1,:)));
antenna = 1;
figure;
for i = 1:6
    temp = eval(var_name(i,:));
    if antenna == 1
        new_data(i,:) =(temp(1,:) + 1i*temp(2,:));
    else
        new_data(i,:) =(temp(3,:) + 1i*temp(4,:));
    end
    
    %         new_data(i,:) =( temp(3,:) + 1i*temp(4,:));
    [find,lags] = xcorr(real(new_data(i,:)),chirp_r);
    len = length(real(new_data(2,:)));
    [value_max,p_max] = max(abs(find(1:3/2*len))); 
    pp(i) = lags(p_max);
    R = IQ_rate/30.72; 
%     res_cor_com(i,:) = xcorr(new_data(i,pp-330000:pp-330000+200000),com_tx_r(1:170000));
    if pp(i) < floor(100000*R)
        pp(i) = pp(i) + len/2;
    end
    res_cor_com(i,:) = xcorr(new_data(i,pp(i)-floor(100000*R):pp(i)+floor((-100000+200000)*R)),com_tx_r(1:floor(200000*R)));
    if(i == 1)
        normal = max(abs(res_cor_com(i,:)));
    end
    plot(abs(res_cor_com(i,:)));
    hold on;
end
legend('1','2','3','4','5','6');
%%
%xcorr data  of synchronized signal
nr = zeros(1,length(new_data(i,pp(i)-floor(100000*R):pp(i)+floor((-100000+200000)*R))));
nr = complex(nr);
for i = 1:6
    nr = nr + new_data(i,pp(i)-floor(100000*R):pp(i)+floor((-100000+200000)*R));
end
res_cor_com(7,:) = xcorr(nr,com_tx_r(1:floor(200000*R)));
normal = max(abs(res_cor_com(7,:)));
res_cor_com(7,:) = (res_cor_com(7,:))./normal;
figure
plot(abs(res_cor_com(7,:)));
% temp_p = int32(length(res_cor_com)/2-10);
% [a,p] = max(abs(res_cor_com(7,temp_p-20:(temp_p+100))));
% la = linspace((1-p)/30.72,(121-p)/30.72,121);
% plot(abs(res_cor_com(7,temp_p-20:(temp_p+100))));
%%
rr_std = zeros(1,length(la));
% delay_std = [0 2.05e-7 2.85e-7 6.6e-7 8.05e-7 9.25e-7];
delay_std_1 = [0 2.0e-7 2.8e-7 6.6e-7 8.0e-7 9.2e-7];
res_std = [1 0.5346 0.7471 0.3718 0.2504 0.1435];
delay = delay_std_1*1e8 + p;
for i = 1:6
    rr_std(delay(i)) = res_std(i);
end
stem(la,rr_std);
axis([-0.5 1.5 0.05 1.05 ]);
legend('Synthesized data with correct','Model standard');
%%

for i = 1:6
    res_new(i,:) = res_cor_com(i,2764801-50:2764830+50);
    
    res_new_interp(i,:) = resample(abs(res_new(i,:)),10000,3072);ds
    res_new_interp_1(i,:) = resample(abs(res_new_interp(i,:)),10,1);
    
    if i == 1
        normal = max(abs(res_new_interp(i,:)));
    end
    res_new_interp(i,:) = abs(res_new_interp(i,:))./normal;
    
    la_interp = linspace((1-165)/100,(424-165)/100,424);
    la_interp_1 = linspace((1-1629)/1e9,(4240-1629)/1e9,4240);
    
    
    %     plot(la_interp, abs(res_new_interp(i,:)));
    %     plot(la_interp,  abs(res_new_interp(i ,:)))
    hold on;
end

%%
figure;
res_new(7,:) = res_cor_com(7,2764801-50:2764830+50);
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
%%
% %normalization h(t)
% for i = 1:6
%     %     ress(i) = cursor_info_3D_tripod_syn(i).Position(2)/cursor_info_3D_tripod_syn(1).Position(2);
%     %     ress(i) = cursor_info_single_101(i).Position(2)/cursor_info_single_101(1).Position(2);
%     %       ress(i) = cursor_info_3D_turn_syn(i).Position(2)/cursor_info_3D_turn_syn(1).Position(2).*(max_s21(1)/max_s21(i));
%     ress(i) = cursor_info_3D_turn_syn(i).Position(2)/cursor_info_3D_turn_syn(1).Position(2);
% end
% res_std = [1 0.5346 0.7471 0.3718 0.2504 0.1435];
% figure;
% stem(ress);
% hold on;
% stem(res_std);