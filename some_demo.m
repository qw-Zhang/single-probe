data = load('corr_data.dat');
Tx = data([1:4],:);
Rx = data([5:8],:);
td = zeros(2,36300);
rd = zeros(2,length(data));
for i = 1:2
    td(i,:) = data(2*i-1,1:36300) + 1j*data(2*i,1:36300);
    rd(i,:) = data(2*i-1+4,:) + 1j*data(2*i+4,:);
end

%%
tx_corr = td(1,:)+td(2,:);
[a,b] = xcorr(rd(1,:),tx_corr);
figure;plot(b,abs(a));
%%
stat = load('stat.dat');

j = 1;
for i = 1:length(stat)
    if stat(i) ~= 2
        res(j) = stat(i);
        j = j + 1;
    end
end
res_stat = sum(res)/length(res);
%%
path(path,'C:\Exprience file\SMART-NI Back Up\CE\ChannelModelMatCodes\SCM');
scmp = scmparset;linkp = linkparset;antp = antparset;
scmp.NumPaths=6;scmp.NumTimeSamples=1;scmp.FixedPdpUsed = 'yes';
[H delay output] = scm(scmp,linkp,antp);
%   H transform to H_4(h11,h12,h21,h22)
H_s = size(H);
if length(H_s) == 3
    H_4 = zeros(4,H_s(3)*1);
    H_4(1,:) = reshape(H(1,1,:),1,H_s(3)*1);
    H_4(2,:) = reshape(H(1,2,:),1,H_s(3)*1);
    H_4(3,:) = reshape(H(2,1,:),1,H_s(3)*1);
    H_4(4,:) = reshape(H(2,2,:),1,H_s(3)*1);
else
    H_4 = zeros(4,H_s(3)*H_s(4));
    H_4_len = length(H_4);
    H_4(1,:) = reshape(H(1,1,:,:),1,H_s(3)*H_s(4));
    H_4(2,:) = reshape(H(1,2,:,:),1,H_s(3)*H_s(4));
    H_4(3,:) = reshape(H(2,1,:,:),1,H_s(3)*H_s(4));
    H_4(4,:) = reshape(H(2,2,:,:),1,H_s(3)*H_s(4));
end
%%
for i = 1:4
    p(i,:)=power(abs(H_4(i,:)),2);
    p(i,:) = p(i,:)/sum(p(i,:));
end

mean(p,1)
%%
tx_temp = load('C:\Exprience file\SMART-NI Back Up\SMART-NI Back Up\Transmitter\Txdata2.dat');
channel_temp = load('C:\Exprience file\SMART-NI Back Up\SMART-NI Back Up\Transmitter\channel_matrix.dat');
%%
channel = zeros(4,length(channel_temp));
for i = 1:4
    channel(i,:) = channel_temp(2*i-1,:) + 1j*channel_temp(2*i,:);
end

tx = zeros(2,length(tx_temp));
rx = zeros(2,length(tx));
for i = 1:2
    tx(i,:) = tx_temp(2*i-1,:) +  1j*tx_temp(2*i,:);
end

rx(1,:) = conv(tx(1,:),channel(1,:),'same');
rx(1,:) = rx(1,:) + conv(tx(2,:),channel(3,:),'same');
rx(2,:) = conv(tx(2,:),channel(4,:),'same');
rx(2,:) = rx(2,:) + conv(tx(1,:),channel(2,:),'same');
%%

%%
figure;hold on;
plot(real(tx(1,:)));
plot(real(rx(1,:)));
%%

scmpar=scmparset;
scmpar.NumTimeSamples=6; 
scmpar.NumPaths=6;
[H tau]=scm(scmpar,linkparset,antparset);
new_H = reshape(H,2,72);
figure;plot(abs(sum(channel,1)))
%%
[a,b] = xcorr(rx(1,:),tx(1,:));
figure;plot(b,abs(a));

%%
path(path,'C:\Exprience file\SMART-NI Back Up\SMART-NI Back Up\Receiver\data_record');
data_4M = load('data_source.dat');
%%
data(1,:) = data_4M(1,:) + 1j*data_4M(2,:);
data(2,:) = data_4M(3,:) + 1j*data_4M(4,:);
%%
for i = 1:1
    data_up(i,:) = resample(data(i,:),25,1);
end
channel  = zeros(1,1000); %resolution 10ns
% channel(1) = 1;
% channel(50) = 0.5;
channel(100) = 0.7; 
channel(500) = 0.5;
channel(700) = 0.2;

data_channel = conv(real(data_up(1,:)),channel);
data_down = resample(data_channel,1,25);

channel_4M = zeros(1,250);
channel_4M(1) = 1;
channel_4M(2) = 0.5;
channel_4M(4) = 0.7;
channel_4M(20) = 0.3;
channel_4M(28) = 0.2;

data_r = conv(data(1,:),channel_4M);
%%
figure;hold on;
plot(real(data_channel));
% plot(real(data_r));
% plot(real(data(1,:)));
%%
[a,b] = xcorr(data_down,data(1,102413:102413+10000));
figure;plot(b,abs(a));
%%
path(path,'C:\Exprience file\SMART-NI Back Up\SMART-NI Back Up\Receiver\data_record');
data_channel = load('data_channel.dat');
data_source = load('data_source.dat');

%%
P = -70:2:-28;
val_ce = zeros(1,length(P));
val = zeros(1,length(P));

BER_pos_ce = -50:2:-28;
BER_val_ce = [1,0.747103,0.530077,0.267538,3.299E-10,3.96499E-13,5.97319E-21,5.16013E-32,7.11731E-47,3.27498E-72,1.54837E-60,3.47163E-95];

BER_pos = -60:2:-28;
BER_val = [1,0.743114,0.0877759,0.0208312,0.06721,3.45456e-32,1.59529e-83,1.28468e-71,5.817e-94,7.36226e-172,1.60695e-128,1.36727e-71,5.8976e-103,5.8976e-103,5.8976e-103,5.8976e-103,5.8976e-103];

for i = 1:length(BER_val_ce)
    val_ce(P == BER_pos_ce(i)) = 1 - BER_val_ce(i);
end
for i = 1:length(BER_val)
    val(P == BER_pos(i)) = 1 - BER_val(i);
end
figure;hold on;
plot(P,val_ce,'Marker','*','LineWidth',2);
plot(P,val,'Marker','o','LineWidth',2);

%%
be = [0,0.00152721,(0.157298+0.328139+0.367368+0.26573)/4,(0.500396+0.498213+0.182822+0.083025)/4,(0.013991+0.175031+0.261477)/3,...
    0.275992,(0.303552+0.465756+0.214556+0.483242+0.517729)/5,(0.619461+0.700947+0.686722+0.655699+0.390936+0.672337)/6,...
    (0.6493+0.792758+0.83075)/3];

ber_1 = [0,0.00152721,0.157298,0.182822,0.261477,0.275992,0.465756,(0.619461+0.700947+0.686722+0.655699+0.390936+0.672337)/6,...
    (0.6493+0.792758+0.83075)/3];
ber_1 = flip(ber_1);
BER_pos_1 = -48:1:-40;
P1 = -70:1:-28;
val_1 = zeros(1,length(P1));
for i = 1:length(ber_1)
    val_1(P1 == BER_pos_1(i)) = 1 - ber_1(i);
    pos_temp = find(P1 == BER_pos_1(i));
end
val_1(pos_temp:length(P1)) = 1;
figure;
plot(P1,val_1,'Marker','*','LineWidth',2);
