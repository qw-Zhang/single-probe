% c = load('D:\study\OTA\expriment_iecas\DATA\DATA0923\chirp_100_30.72_80.dat');
% c = load('D:\study\OTA\expriment_iecas\DATA\DATA0923\chirp_100_90dbm.dat');
c_s = load('D:\study\OTA\expriment_iecas\DATA\DATA0923\chirp_std.dat');
c_std = c_s + 1i*c_s;
c_20 = load('chirp_100_20_80_1.dat');
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
signal_ex_len = length(c)/(2*rate); %I(0) Q(0) I(1) Q(1) data formal of original .dat file
signal_piont_ex = zeros(signal_ex_len,rate);

for i = 1:length(c)/2
    c_r(i) = c(2*i-1) + 1i*c(2*i);  %extra the oversampling data from original .dat file
end

%if we should "delayseq" function, the formal data should array as one row
%for a group data
%loop for rate
for r = 1:rate
    j = 1;
    %extra same sampling rate data from oversampling data
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
%avearage rate times statisitcal result every row
%mean's second parameter default is 1->row,2->column
p_mean = mean(p);
figure;hist(p_mean);

%%
%these function for .dat which too big. the data perform 2x2200000 of .dat
%file
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
%%
r_p_mean = p_mean - 1.001645e4;
histfit(r_p_mean,5,'normal');   %firstly generate the fitting line
h = findobj(gca,'Type','line'); %get the fitting line from histfitgram
xdata=get(h,'XData');
ydata=get(h,'YData') ;
yy = ydata/100;     %adjust the data of fitting line
figure;bar(ht, hf/sum(hf));hold on;plot(xdata,yy);
%%
% piont_move = 17;
piont_move = mean(p,2);
for r = 1:rate
    for i = 1:signal_ex_len
        signal_piont_ex(j,r) = B_r(rate*(i-1)+r);
        j = j + 1;
    end
    signal_piont_ex_in(:,r) = resample(signal_piont_ex(:,r),rate,1);
    signal_piont_ex_in_m(:,r) = delayseq(signal_piont_ex_in(:,r),piont_move(r)-10000);
end

for r = 1:rate
    j = 1;
    for i = 1:signal_ex_len
        signal_piont_ex_m(j,r) = signal_piont_ex_in_m(rate*(i-1)+1,r);
        j = j + 1;
    end
    for i = 1:101
        temp_m(:,i) = signal_piont_ex_m((i-1)*10000 + 1:(i-1)*10000 + 10000,r);
        [a,b] = xcorr(temp_m(:,i),c_s(1:5000));
        [v_temp_m(2*r,i),v_temp_m(2*r-1,i)] = max(abs(a));
        p_m(r,i) = v_temp_m(2*r-1,i);
    end
    
end

p_mean_m = mean(p_m);
figure;histogram(p_mean_m,20);