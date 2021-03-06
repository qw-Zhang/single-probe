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
% dir_dic = 'D:\study\OTA\expriment_iecas\DATA\DATA0811\3D_turn';

for i = 1:6
    temp_ch = num2str(i);
%     name(i,:) = [dir_dic,'\data',temp_ch,'_chirp',temp_ch_pow,'dbm','.dat'];
    name(i,:) = [dir_dic,'\data',temp_ch,'_c',temp_ch_pow,'dbm','.dat'];
%     var_name(i,:) = ['data',temp_ch,'_chirp_10dbm'];
    var_name(i,:) = ['data',temp_ch,'_c10dbm'];
    load(name(i,:));
end

f0 = 10;f1 = 2000;
t = linspace(1,10000,10000);
chirp = sin(2 * pi*(f0*t / 10000 + (f1 / 2)*power(t / 10000, 2)));
%%
% %gen VNA calibrate coefficient
% power_map = [-15 -10 0 5 10];
% figure;
% for i = 1:6
%     for j = 1:3
%         name_p = [dir_dic,'\VNA_',num2str(angle_map(i)),'_90_',num2str(power_map(j)),'.dat'];
%         p_vna(:,(5*j-4):(5*j)) = load(name_p);
%         s21(:,j) = p_vna(:,(5*j-4+1)) + 1i*p_vna(:,(5*j-4+2));
%         s31(:,j) = p_vna(:,(5*j-4+3)) + 1i*p_vna(:,(5*j-4+4));
%         ifft_s21(:,j) = ifft(s21(:,j));
%         ifft_s31(:,j) = ifft(s31(:,j));
%     end
%     mean_ifft_s21(:,i) = mean(ifft_s21,2);
%     mean_ifft_s31(:,i) = mean(ifft_s31,2);
%     max_s21(i) = max(abs(mean_ifft_s21(:,i)));
%     max_s31(i) = max(abs(mean_ifft_s31(:,i)));
%     plot(abs(mean_ifft_s31(:,i)));
%     hold on;
% end
% legend('1','2','3','4','5','6');

%%
%gen chirp calibrate cofficient
R = 1;
antenna = 1;correct = 0;
new_data_chirp = zeros(6,length(eval(var_name(1,:))));
len = length(new_data_chirp);
for i = 1:6
    temp = eval(var_name(i,:));
    if antenna == 1
        new_data_chirp(i,:) =(temp(1,:) + 1i*temp(2,:));
    else
        new_data_chirp(i,:) =(temp(3,:) + 1i*temp(4,:));
    end
    [find,lags] = xcorr((new_data_chirp(i,:)),chirp);
%     figure;    plot(lags,find);
    [value_max(i),p_max(i)] = max(abs(find(1:3/2*len)));
	pp_value_max(i) = lags(p_max(i));
    if pp_value_max(i) < floor(100000*R)
        pp_value_max(i) = pp_value_max(i) + 1230800;
        value_max(i) = abs(find(pp_value_max(i)+2764800));
    end
end
%%
%xcorr data  of respective angle
new_data = zeros(7,length(eval(var_name(1,:))));
r = zeros(7,2*length(eval(var_name(1,:)))-1);
com_tx = Txdata(1,:) + 1i*Txdata(2,:);
% w = hamming(length(data110dbm(1,:)));

figure;
for i = 1:6
    temp = eval(var_name(i,:));
    if antenna == 1
        if correct == 1
            new_data(i,:) =(value_max(1)/value_max(i)).*( temp(1,:) + 1i*temp(2,:));
        else
            new_data(i,:) =(temp(1,:) + 1i*temp(2,:));
        end
    else
        if correct == 1
            new_data(i,:) =(value_max(1)/value_max(i)).*( temp(3,:) + 1i*temp(4,:));
        else
            new_data(i,:) =(temp(3,:) + 1i*temp(4,:));
        end
    end
    
    %         new_data(i,:) =( temp(3,:) + 1i*temp(4,:));
%     [res_cor_com(i,:),lags] = xcorr(new_data(i,:),com_tx);
    if pp_value_max(i) < floor(100000*R)
        pp_value_max(i) = pp_value_max(i) + len/2;
    end
    res_cor_com(i,:) = xcorr(new_data(i,pp_value_max(i)-floor(100000*R):pp_value_max(i)+floor((-100000+200000)*R)),com_tx(1:floor(200000*R)));
    %     [res_cor_com(i,:),lags] = xcorr(new_data(i,:),com_tx);
    if(i == 1)
        normal = max(abs(res_cor_com(i,:)));
    end
    plot(abs(res_cor_com(i,:)));
    hold on;
end
legend('1','2','3','4','5','6');
%%
%xcorr data  of synchronized signal
nr = zeros(1,length(new_data(i,pp_value_max(i)-floor(100000*R):pp_value_max(i)+floor((-100000+200000)*R))));
nr = complex(nr);
nr_res = zeros(1,length(res_cor_com(1,:)));
for i = 1:6
    nr = nr + new_data(i,pp_value_max(i)-floor(100000*R):pp_value_max(i)+floor((-100000+200000)*R));
    nr_res = nr_res + res_cor_com(i,:);
end
res_cor_com(7,:) = xcorr(nr,com_tx(1:floor(200000*R)));
normal = max(abs(res_cor_com(7,:)));
res_cor_com(7,:) = (res_cor_com(7,:))./normal;
figure
plot(abs(res_cor_com(7,:)));
mm = max(abs(nr_res));
% nr_res = nr_res./mm;
% figure;plot(abs(nr_res));

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
[start_value,start_index] = max(res_cor_com(1,:));
res_new = zeros(7,length(start_index-50:start_index+50));
for i = 1:6
    res_new(i,:) = res_cor_com(i,start_index-50:start_index+50);
    
    res_new_interp(i,:) = resample(abs(res_new(i,:)),10000,3072);
    res_new_interp_1(i,:) = resample(abs(res_new_interp(i,:)),10,1);
    
    if i == 1
        normal = max(abs(res_new_interp(i,:)));
    end
    res_new_interp(i,:) = abs(res_new_interp(i,:))./normal;
    
    la_interp = linspace((1-165)/100,(329-165)/100,329);
%     la_interp_1 = linspace((1-1629)/1e9,(4240-1629)/1e9,4240);
    
    
    %     plot(la_interp, abs(res_new_interp(i,:)));
    %     plot(la_interp,  abs(res_new_interp(i ,:))) 
%     hold on;
end

%%
figure;
res_new(7,:) = res_cor_com(7,start_index-50:start_index+50);
res_new_interp(7,:) = resample(abs(res_new(7,:)),10000,3072);
res_new_interp(7,:) = res_new_interp(7,:)./max(res_new_interp(7,:));
res_new_interp_1(7,:) = resample(abs(res_new_interp(7,:)),10,1);
plot(la_interp,  abs(res_new_interp(7 ,:)),'LineWidth',5);hold on;
xlabel('Delay [microsecond]');ylabel('Normalization Amplitude')
%%
rr_std = zeros(1,length(la_interp));
% delay_std = [0 2.05e-7 2.85e-7 6.6e-7 8.05e-7 9.25e-7];
delay_std_1 = [0 2e-7 2.8e-7 6.6e-7 8.0e-7 9.2e-7];
res_std = [1 0.5346 0.7471 0.3718 0.2504 0.1435];
delay = delay_std_1*1e8 + 165;
for i = 1:6
    rr_std(delay(i)) = res_std(i);
end
stem(la_interp,rr_std);
axis([-0.5 1.5 0.05 1.05 ]);
legend('Synthesized data with correct','Model standard');