%%
com1 = zeros(4,10,1000);
for n = 1:4
    str = num2str(n);
    dir = 'D:\study\OTA\expriment_ucas\0309\data';
    filename = [dir,str,'.txt'];
    data = load(filename);
    
    batch_size = 1000;
    it_size = length(data) / batch_size;
    amp_mean = zeros(1,it_size);
    ang_mean = zeros(1,it_size);
    data_size = size(data);
    
    for i = 1:it_size
        start = (i - 1) * batch_size + 1;
        stop = i * batch_size;
        com1(n,i,:) = data( start:stop ,1) + 1i * data( start:stop ,2);
        amp = abs(squeeze(com1(n,i,:)));
        ang = angle(squeeze(com1(n,i,:)));
        amp_mean(i) = mean(mag2db(amp));
        ang_mean(i) = mean(ang);
        %     [mean(mag2db(amp)), mean(ang)]
    end
    figure;
    plot(ang_mean);
    
end
%%
data_pro = zeros(8,9,1000);
loc = 90.525:-8:26.525;
data_pro(1:4,:,:) = com1(1:4,2:10,:);
%%
n = 8;
str = num2str(n);
dir = 'D:\study\OTA\expriment_ucas\0309\data';
filename = [dir,str,'.txt'];
data = load(filename);

batch_size = 1000;
it_size = length(data) / batch_size;
amp_mean = zeros(1,it_size);
ang_mean = zeros(1,it_size);
data_size = size(data);

for i = 1:it_size
    start = (i - 1) * batch_size + 1;
    stop = i * batch_size;
    com1(i,:) = data( start:stop ,1) + 1i * data( start:stop ,2);
    amp(i,:) = abs(squeeze(com1(i,:)));
    ang(i,:) = angle(squeeze(com1(i,:)));
    amp_mean(i) = mean(mag2db(amp(i,:)));
    ang_mean(i) = mean(ang(i,:));
    %     [mean(mag2db(amp)), mean(ang)]
end
data_pro(n,:,:) = com1;
% figure;
% plot(ang_mean);

%%
% data_pro 8x9x1000 angle x location x samples
% load('D:\study\OTA\expriment_ucas\0309\data_fre_res.mat');
% loc = 90.525:-8:26.525; % slide resolution

% data_pro 8x13x1000 angle x location x samples
load('D:\study\OTA\expriment_ucas\0310\data_pro_micro.mat');
loc = 98.525:-8:2.525; % slide resolution

data_size = size(data_pro);
dim2 = data_size(2);
% x axis
fc = 2.45e9;
speed_of_light=2.99792458e8;
lambda = speed_of_light / fc;
for i = 1: (floor(dim2 / 2) + 1)
    dis(i) = ( loc(i) - loc(dim2 + 1 - i) ) * 0.125 * 1e-2;
end
dis_lam = dis ./ lambda;
dis_lam = fliplr(dis_lam);

% m1
for i = 1:data_size(1)
    
    for j = 1:int32(dim2/2)
        ant1 = data_pro(i, j ,:);
        ant2 = data_pro(i, (dim2 + 1 - j) ,:);
        
        ang_ant1 = angle(ant1);
        ang_ant2 = angle(ant2);
        ang_diff(i,j) = exp(1i*mean(ang_ant2 - ang_ant1));
        
        corr_temp = corrcoef(ant1,ant2);
        corr(i,j) = corr_temp(2);
        mean_m1_temp(i,j) = mean(ant1 .* conj(ant2));
    end
end
mean_m1 = sum(mean_m1_temp);
mean_m1 = fliplr(mean_m1);
mean_m1 = mean_m1 ./ mean_m1(1);

corr_ang = sum(ang_diff,2);

figure;
plot(dis_lam,  abs(mean_m1));
%%
m2_data = squeeze(sum(data_pro,1));
for j = 1:int32(dim2/2)
    ant1_m2(j,:) = m2_data(j ,:);
    ant2_m2(j,:) = m2_data((dim2 + 1 - j) ,:);
    corr_temp = corrcoef(ant1_m2(j,:),ant2_m2(j,:));
    corr_m2(j) = corr_temp(2);
    mean_m2(j) = mean(ant1_m2(j,:) .* conj(ant2_m2(j,:)) );
end

mean_m2 = fliplr(mean_m2);
mean_m2 = mean_m2 ./ mean_m2(1);

corr_m2 = fliplr(corr_m2);
% figure;
hold on;
plot(dis_lam,  abs(mean_m2));
%%
corr1 = fliplr(corr);
c = sum(corr1,1);
cc = c ./ c(1);
figure;
plot(dis_lam,  abs(cc));
%%
% 0310 micro + macro
com1 = zeros(8,13,1000);
for n = 1:8
    str = num2str(n);
    dir = 'D:\study\OTA\expriment_ucas\0310\data';
    filename = [dir,str,'_a.txt'];
    data = load(filename);
    
    batch_size = 1000;
    it_size = length(data) / batch_size;
    amp_mean = zeros(1,it_size);
    ang_mean = zeros(1,it_size);
    data_size = size(data);
    
    for i = 1:it_size
        start = (i - 1) * batch_size + 1;
        stop = i * batch_size;
        com1(n,i,:) = data( start:stop ,1) + 1i * data( start:stop ,2);
        amp = abs(squeeze(com1(n,i,:)));
        ang = angle(squeeze(com1(n,i,:)));
        amp_mean(i) = mean(mag2db(amp));
        ang_mean(i) = mean(ang);
        %     [mean(mag2db(amp)), mean(ang)]
    end
    
    noun = 0;
%     figure;
%     plot(ang_mean);
%     
end
%%
% deal the data5_i...
data_pro_0310 = zeros(8,13,1000);
data_pro_0310(1:4,:,:) = com1(1:4,1:13,:);
data_pro_0310(5,1:4,:) = com1(5,1:4,:);
data_pro_0310(5,5:13,:) = com1(5,6:14,:);
data_pro_0310(6:8,:,:) = com1(6:8,1:13,:);
%%
data_pro_0310_macro = com1;