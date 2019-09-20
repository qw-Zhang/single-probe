clear;
power_map = [-15 -10 0 5 10];
angle_map = [1 -13 146 -30 -11 -1];
 figure;
for i = 1:3
    for j = 1:5
%         name_p = ['VNA_',num2str(angle_map(i)),'_90_',num2str(power_map(j)),'.dat'];
        name_p = ['p',num2str(i),'_',num2str(power_map(j)),'dbm.dat'];
%         name_p = ['D:\study\OTA\expriment_iecas\DATA\DATA0809\turn\p',num2str(i),'_',num2str(power_map(j)),'dbm.dat'];
        p_vna(:,(5*j-4):(5*j)) = load(name_p);
        s21(:,j) = p_vna(:,(5*j-4+1)) + 1i*p_vna(:,(5*j-4+2));
        s31(:,j) = p_vna(:,(5*j-4+3)) + 1i*p_vna(:,(5*j-4+4));
        ifft_s21(:,j) = ifft(s21(:,j));
        ifft_s31(:,j) = ifft(s31(:,j));
    end
     mean_ifft_s21(:,i) = mean(ifft_s21,2);
     mean_ifft_s31(:,i) = mean(ifft_s31,2);
     plot(abs(mean_ifft_s21(:,i)));
     hold on;
end

for i = 1:6
    for j = 1:5
%         name_p = ['VNA_',num2str(angle_map(i)),'_90_',num2str(power_map(j)),'.dat'];
%         name_p = ['D:\study\OTA\expriment_iecas\DATA\DATA0809\tripod\p',num2str(i),'_',num2str(power_map(j)),'dbm.dat'];
        name_p = ['D:\study\OTA\expriment_iecas\DATA\DATA0809\turn\p',num2str(i),'_',num2str(power_map(j)),'dbm.dat'];
        p_vna(:,(5*j-4):(5*j)) = load(name_p);
        s21(:,j) = p_vna(:,(5*j-4+1)) + 1i*p_vna(:,(5*j-4+2));
        s31(:,j) = p_vna(:,(5*j-4+3)) + 1i*p_vna(:,(5*j-4+4));
        ifft_s21(:,j) = ifft(s21(:,j));
        ifft_s31(:,j) = ifft(s31(:,j));
    end
     mean_ifft_s21(:,i) = mean(ifft_s21,2);
     mean_ifft_s31(:,i) = mean(ifft_s31,2);
     plot(abs(mean_ifft_s21(:,i)));
     hold on;
end
 legend('t1','t2','t3','1','2','3','4','5','6');
% plot(abs(ifft_s21));
% hold on;
% plot(abs(ifft_s31));