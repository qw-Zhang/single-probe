name = 'D:\study\OTA\expriment_iecas\DATA\20190810\multi_tripod\1-15dbm.dat';
    var_name = 'X1_15dbm';
    load(name);
    
    %%
%xcorr data  of respective angle
new_data_re = zeros(6,2764800);
r = zeros(7,5529599);
com_tx = Txdata(1,:) + 1i*Txdata(2,:);
% w = hamming(length(data110dbm(1,:)));
figure;
for i = 1:1
    temp = eval(var_name(i,:));
    new_data(i,:) = temp(1,:) + 1i*temp(2,:);
    res_cor_com(i,:) = xcorr(new_data(i,:),com_tx);
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
    plot(abs(res_cor_com(i,:)));
    hold on;
end
legend('1','2','3','4','5','6');
    %%
power_map = [-15 -10 0 5 10];
figure;
for i = 1:3
    for j = 1:5
        name_p = ['D:\study\OTA\expriment_iecas\DATA\20190810\p',num2str(i),'_',num2str(power_map(j)),'dbm.dat'];
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
%normalization h(t)
for i = 1:6
%     ress(i) = cursor_info_multi_110(i).Position(2)/cursor_info_multi_110(1).Position(2)*(amp_correct_multi(1)/amp_correct_multi(i));
    ress(i) = cursor_info_multi_10(i).Position(2)/cursor_info_multi_10(1).Position(2);
%     ress(i) = cursor_info(i).Position(2)/cursor_info(1).Position(2);
end
res_std = [1 0.5346 0.7471 0.3718 0.2504 0.1435];
figure;
stem(ress);
hold on;