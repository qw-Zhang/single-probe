power_map = [-15 -10 0 5 10];
angle_map = [1 -13 146 -30 -11 -1];
figure;
for i = 1:6
   
%         name_p = ['VNA_',num2str(angle_map(i)),'_90_',num2str(power_map(j)),'.dat'];
        name_p = ['p',num2str(i),'_',num2str(10),'dbm.dat'];
        p_vna = load(name_p);
        s21 = p_vna(:,(1)) + 1i*p_vna(:,(2));
        s31 = p_vna(:,(3)) + 1i*p_vna(:,(4));
        ifft_s21 = ifft(s21);
        ifft_s31 = ifft(s31);
   
%      mean_ifft_s21(:,i) = mean(ifft_s21,2);
%      mean_ifft_s31(:,i) = mean(ifft_s31,2);
     plot(abs(ifft_s31));
     hold on;
end
legend('1','2','3','4','5','6');
% plot(abs(ifft_s21));
% hold on;
% plot(abs(ifft_s31));