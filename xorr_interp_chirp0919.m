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
%     name(i,:) = [dir_dic,'\data',temp_ch,'_chirp',temp_ch_pow,'dbm','.dat'];
    name(i,:) = [dir_dic,'\data',temp_ch,'_c',temp_ch_pow,'dbm','.dat'];
%     var_name(i,:) = ['data',temp_ch,'_chirp_10dbm'];
    var_name(i,:) = ['data',temp_ch,'_c10dbm'];
    load(name(i,:));
end

f0 = 10;f1 = 2000;
t = linspace(0,9999,10000);
chirp = sin(2 * pi*(f0*t / 10000 + (f1 / 2)*power(t / 10000, 2)));

 %%
% %using interpolation 
% R = 1;
% antenna = 2;correct = 1;
% new_data_chirp = zeros(6,length(eval(var_name(1,:))));
% len = length(new_data_chirp);
% rate = 10; %the interpolation rate
% 
% s_ex = zeros(6,2764800);
% 
% for m = 1:1
% %      k = -rate/2;
% %     k = 1;
%     temp = eval(var_name(m,:));
%     if antenna == 1
%         new_data_chirp(m,:) =(temp(1,:) + 1i*temp(2,:));
%     else
%         new_data_chirp(m,:) =(temp(3,:) + 1i*temp(4,:));
%     end
%     s_interp = resample(new_data_chirp(m,:),rate,1);
% %     s_interp = interp1(new_data_chirp(m,:),rate);
%     s_ex_len = length(s_interp)/rate;
%     
%     for n = 1:rate
%         j = 1;
%         for i = 1:s_ex_len
%             s_ex(m,j) = s_interp(rate*(i-1)+n);
%             j = j + 1;
%         end
%         [a,b] = xcorr(s_ex(m,:),chirp);
%         [p(2,n),p(1,n)] = max(abs(a(100000:3/2*len)));
%         
% %         p(1,n) = p(1,n) - length(s_ex);
%         p(1,n) = b(p(1,n));
%         %     figure;plot(b,a);    
%         if p(1,n) < floor(100000*R)
%             p(1,n) = p(1,n) + 1229800; 
%             p(2,n) = abs(a(p(1,n)+2763800));
%         end
% %         k = k + 1;
%     end
%     
%     [value_index,index] = max(p(2,:));
%     i = 1;
%     temp_max = p(2,1);index_max = p(1,1);
%     while (i <= rate) && (p(1,1) == p(1,i))
%         if temp_max < p(2,i)
%             temp_max = p(2,i);
%             index_max = p(1,i);
%         end
%         i = i + 1;
%     end
% %     rr(1,m) = p(1,index);
%     rr(1,m) = index_max;
% %     rr(2,m) = value_index;
%     rr(2,m) = temp_max;
%     
% end
% pp_value_max = rr(1,:);
% value_max = rr(2,:);
%%
%using interpolation
R = 1;
antenna = 2;correct = 1;
new_data_chirp = zeros(length(eval(var_name(1,:))),6);
len = length(new_data_chirp);
rate = 20; %the interpolation rate

s_ex = zeros(2764800,6);

for m = 1:6
    temp = eval(var_name(m,:));
    if antenna == 1
        new_data_chirp(:,m) =(temp(1,:) + 1i*temp(2,:));
    else
        new_data_chirp(:,m) =(temp(3,:) + 1i*temp(4,:));
    end
    
    s_interp(:,m) = resample(new_data_chirp(:,m),rate,1);
    s_ex_len = length(s_interp(:,1))/rate;
    
    for n = 1:rate
        j = 1;
        for i = 1:s_ex_len
            s_ex(j,m) = s_interp(rate*(i-1)+n,m);
            j = j + 1;
        end
        [a,b] = xcorr(s_ex(:,m),chirp);
        [p(2*m,n),p(2*m-1,n)] = max(abs((a((s_ex_len + 100000):3/2*len))));
%         p(2*m-1,n) = b(p(2*m-1,n));
    end
    
end

s_result_ex = zeros(2764800,6);
for m = 1:6
    [v_max(2,m),v_max(1,m)] = max(p(2*m,:));
    p_move = rate + 2 - v_max(1,m);
    s_interp_p = delayseq(s_interp(:,m),p_move);
    j = 1;
    for i = 1:s_ex_len
        s_result_ex(j,m) = s_interp_p(rate*(i-1) + 1);
        j = j + 1;
    end
    
end
%%
%gen chirp calibrate cofficient
% new_data_chirp = s_result_ex;
for i = 1:6
    [find,lags] = xcorr((s_result_ex(:,i)),chirp);
%     figure;    plot(lags,find);
    [value_max(i),p_max(i)] = max(abs(find(1:3/2*len)));
	pp_value_max(i) = lags(p_max(i));
    if pp_value_max(i) < floor(100000*R)
        pp_value_max(i) = pp_value_max(i) + 1230800;
        value_max(i) = abs(find(pp_value_max(i)+2764800));
    end
end
% for m = 1:6
% pp_value_max(m) = p(2*m-1,v_max(1,m));
% value_max(m) = p(2*m,v_max(1,m));
% end

%%
%xcorr data  of respective angle
% new_data = zeros(7,length(eval(var_name(1,:))));
signal_p = s_result_ex;
new_data = zeros(length(signal_p),7);
r = zeros(2*length(eval(var_name(1,:)))-1,7);
com_tx = Txdata(1,:) + 1i*Txdata(2,:);

figure;
for i = 1:6
%     temp = eval(var_name(i,:));
    if antenna == 1
        if correct == 1
            new_data(:,i) =(value_max(1)/value_max(i)).*(signal_p(:,i));
        else
            new_data(:,i) =(signal_p(:,i));
        end
    else
        if correct == 1
            new_data(:,i) =(value_max(1)/value_max(i)).*(signal_p(:,i));
        else
            new_data(:,i) =(signal_p(:,i));
        end
    end
    
%     if pp_value_max(i) < floor(100000*R)
%         pp_value_max(i) = pp_value_max(i) + len/2;
%     end
%     res_cor_com(i,:) = xcorr(new_data(i,pp_value_max(i)-floor(100000*R):pp_value_max(i)+floor((-100000+200000)*R)),com_tx(1:floor(200000*R)));
    res_cor_com(i,:) = xcorr(new_data(pp_value_max(i):pp_value_max(i)+floor((-100000+200000)*R),i),com_tx(1:floor(200000*R)));
    if(i == 1)
        normal = max(abs(res_cor_com(i,:)));
    end
    plot(abs(res_cor_com(i,:)));
    hold on;
end
legend('1','2','3','4','5','6');
%%
%xcorr data  of synchronized signal
nr = zeros(length(new_data(pp_value_max(i):pp_value_max(i)+floor((-100000+200000)*R))),1);
nr = complex(nr);
nr_res = zeros(1,length(res_cor_com(1,:)));
for i = 1:6
    nr = nr + new_data(pp_value_max(i):pp_value_max(i)+floor((-100000+200000)*R),i);
    nr_res = nr_res + res_cor_com(i,:);
end
res_cor_com(7,:) = xcorr(nr,com_tx(1:floor(200000*R)));
normal = max(abs((res_cor_com(7,:))));
res_cor_com(7,:) = (res_cor_com(7,:))./normal;
figure
plot(abs((res_cor_com(7,:))));
mm = max(abs(nr_res));
nr_res = nr_res./mm;
% figure;
hold on;
plot(abs(nr_res));

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
%     hold on;
end

%%
% figure;
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