%%
clear;
%%
%loading data
clear name var_name;
Txdata = load('D:\study\OTA\expriment_iecas\DATA\DATA0730\Tx_data\Txdata.dat');
temp_ch_pow = '10';
angle_map = [1 -13 146 -30 -11 -1];
% dir_dic = 'D:\study\OTA\expriment_iecas\DATA\DATA0810\multi_tripod';
% dir_dic = 'D:\study\OTA\expriment_iecas\DATA\DATA0810\tripod';
% dir_dic = 'D:\study\OTA\expriment_iecas\DATA\DATA0811\2D_turn';
% dir_dic = 'D:\study\OTA\expriment_iecas\DATA\DATA0811\3D_tripod';
% dir_dic = 'D:\study\OTA\expriment_iecas\DATA\DATA0811\3D_turn';
dir_dic = 'D:\study\OTA\expriment_iecas\DATA\DATA0824\3D_1';
for i = 1:6
    temp_ch = num2str(i);
    name(i,:) = [dir_dic,'\data',temp_ch,temp_ch_pow,'dbm2','.dat'];
    var_name(i,:) = ['data',temp_ch,'10dbm2'];
    load(name(i,:));
end

%%
power_map = [-15 -10 0 5 10];
figure;
for i = 1:6
    for j = 1:3
        name_p = [dir_dic,'\VNA_',num2str(angle_map(i)),'_90_',num2str(power_map(j)),'.dat'];
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
legend('1','2','3','4','5','6');

%%
%xcorr data  of respective angle
new_data_re = zeros(6,length(eval(var_name(1,:))));
r = zeros(7,2*length(eval(var_name(1,:)))-1);
com_tx = Txdata(1,:) + 1i*Txdata(2,:);
% w = hamming(length(data110dbm(1,:)));
antenna = 2;correct = 0;
figure;
for i = 1:6
    temp = eval(var_name(i,:));
    if antenna == 1
        if correct == 1
            new_data(i,:) =(max_s31(1)/max_s31(i)).*( temp(1,:) + 1i*temp(2,:));
        else
            new_data(i,:) =(temp(1,:) + 1i*temp(2,:));
        end
    else
        if correct == 1
            new_data(i,:) =(max_s21(1)/max_s21(i)).*( temp(3,:) + 1i*temp(4,:));
        else
            new_data_1(i,:) =(temp(3,:) + 1i*temp(4,:));
        end
    end
    
    %         new_data(i,:) =( temp(3,:) + 1i*temp(4,:));
    [res_cor_com_1(i,:),lags] = xcorr(new_data_1(i,:),com_tx);
    %     [res_cor_com(i,:),lags] = xcorr(new_data(i,:),com_tx);
    if(i == 1)
        normal = max(abs(res_cor_com_1(i,:)));
    end
    %     res_cor_com(i,:) = abs(res_cor_com(i,:))./normal;
%     res_cor_com(i,:) = abs(res_cor_com(i,:));
%     res_cor_com_interp(i,:) = resample(abs(res_cor_com(i,:)),10000,3072);
%     lags_interp = resample(lags,10000,3072);
%     temp_p = int32(length(res_cor_com_interp)/2-10);
    %     plot(r1(i,temp_p:(temp_p+100)));
%         plot(abs(res_cor_com_interp(i,temp_p-30:temp_p+200)));
%     lag_interp = linspace(lags(1),1799997*100-6,17999997);
    plot(abs(res_cor_com_1(i,:)));
    hold on;
end
legend('1','2','3','4','5','6');
%%
la_interp = linspace((1-162)/100,(424-162)/100,424);
figure;
for i=1:6
    new(i,:) = res_cor_com_1(i,2764801-50:2764830+50);
    new_interp(i,:) = resample(abs(new(i,:)),10000,3072);
    if i == 1
        nn = max(new_interp(i,:));
    end
    new_interp(i,:) = new_interp(i,:)./nn;
    plot(la_interp,  abs(new_interp(i ,:)),'LineWidth',5);
    hold on;
end
xlabel('Delay [microsecond]');ylabel('Normalization Amplitude')
%%
%xcorr data  of synchronized signal
nr = zeros(1,length(temp(1,:)));
nr = complex(nr);
for i = 1:6
    nr = nr + new_data(i,:);
end
res_cor_com_1(7,:) = xcorr(nr,com_tx);
normal = max(abs(res_cor_com_1(7,:)));
res_cor_com_1(7,:) = abs(res_cor_com_1(7,:))./normal;
figure
temp_p = int32(length(res_cor_com_1)/2-10);
[a,p] = max(abs(res_cor_com_1(7,temp_p-20:(temp_p+100))));
la = linspace((1-p)/30.72,(121-p)/30.72,121);
plot(abs(res_cor_com_1(7,temp_p-20:(temp_p+100))));
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
    res_new(i,:) = res_cor_com_1(i,2764801-50:2764830+50);
    
    res_new_interp(i,:) = resample(abs(res_new(i,:)),10000,3072);
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
res_new(7,:) = res_cor_com_1(7,2764801-50:2764830+50);
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