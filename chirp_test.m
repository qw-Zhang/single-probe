% c = load('D:\study\OTA\expriment_iecas\DATA\DATA0923\chirp_100_30.72_80.dat');
% c = load('D:\study\OTA\expriment_iecas\DATA\DATA0923\chirp_100_90dbm.dat');
c_s = load('D:\study\OTA\expriment_iecas\DATA\DATA0923\chirp_std.dat');
c_std = c_s + 1i*c_s;
c_20 = load('chirp_100_20_90_1.dat');
B=c_20;B(:,1)=[];B;
%%
% % cc = resample(c_std,10,1);
% c_r = zeros(1,length(c)/2);
% for i = 1:length(c)/2
%     c_r(i) = c(2*i-1) + 1i*c(2*i);
% end
% 
% for i = 1:101
%     temp(i,:) = c_r((i-1)*10000 + 1:(i-1)*10000 + 10000);
%     [a,b] = xcorr(temp(i,:),c_s(1:5000));
%     [v_temp(2,i),v_temp(1,i)] = max(abs(a));
% end

%%
rate = 10;
signal_ex_len = length(c)/(2*rate);
signal_piont_ex = zeros(signal_ex_len,rate);

for i = 1:length(c)/2
    c_r(i) = c(2*i-1) + 1i*c(2*i);
end

for r = 1:rate
    j = 1;
    for i = 1:signal_ex_len
        signal_piont_ex(j,r) = c_r(rate*(i-1)+r);
        j = j + 1;
    end
    for i = 1:101
        temp(:,i) = signal_piont_ex((i-1)*10000 + 1:(i-1)*10000 + 10000,r);
        [a,b] = xcorr(temp(:,i),c_s(1:5000));
        [v_temp(2*r,i),v_temp(2*r-1,i)] = max(abs(a));
        p(r,i) = v_temp(2*r-1,i);
    end
    
end

p_mean = mean(p);
figure;hist(p_mean);

%%
rate = 20;
signal_ex_len = length(B)/rate;
signal_piont_ex = zeros(signal_ex_len,rate);

for i = 1:length(B)
    B_r(i) = B(1,i) + 1i*B(2,i);
end

for r = 1:rate
    j = 1;
    for i = 1:signal_ex_len
        signal_piont_ex(j,r) = B_r(rate*(i-1)+r);
        j = j + 1;
    end
    for i = 1:101
        temp(:,i) = signal_piont_ex((i-1)*10000 + 1:(i-1)*10000 + 10000,r);
        [a,b] = xcorr(temp(:,i),c_s(1:5000));
        [v_temp(2*r,i),v_temp(2*r-1,i)] = max(abs(a));
        p(r,i) = v_temp(2*r-1,i);
    end
    
end

p_mean = mean(p);
figure;histogram(p_mean,50);