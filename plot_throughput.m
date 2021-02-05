load('BER_data.mat');

P = -64:2:-36;
val_ce = zeros(1,length(P));
val = zeros(1,length(P));
BER_pos = -60:2:-28;
for i = 1:length(BER_val)
    val(P == BER_pos(i)) = 1 - BER_val(i);
end

BER_pos_ce = -48:1:-36;
P_1 = -65:1:-35;
val_1 = zeros(1,length(P_1));
for i = 1:length(BER_ce)
    val_1(P_1 == BER_pos_ce(i)) = 1 - BER_ce(i);    
    pos_temp = find(P_1 == BER_pos_ce(i));
end
val_1(pos_temp:length(P_1)) = 1;

P_new = -72:0.1:-45;
BER_theory = berawgn(P_new+72, 'psk', 4,'nondiff');
BER_theory_ce = berfading(P_new+50, 'psk', 4, 1, 500);


val_sp = zeros(1,length(P_new));
val_ce_sp = val_sp;
for i = 1:length(P_new)
    val_sp(i) = power(1 - BER_theory(i),400);   
    val_ce_sp(i) = power(1 - BER_theory_ce(i),100);
end
% val_sp = pchip(P,val,P_new);
% val_ce_sp = pchip(P_1,val_1,P_new);

figure;hold on;
plot(P-1.5,val,'o','MarkerSize',8,'LineWidth',3);
% plot(-60:2:-47,th_measure,'*','MarkerSize',8,'LineWidth',3);
plot(P_new+9,val_sp,'LineWidth',3);
% plot(P_1,val_1*100,'Marker','*','LineWidth',2);
plot(P_new,val_ce_sp*100,'LineWidth',3);
%%
P = -64:2:-36;
val_ce = ones(1,length(P));
val = ones(1,length(P));
BER_pos = -60:2:-28;
for i = 1:length(BER_val)
    val(P == BER_pos(i)) = BER_val(i);
end

BER_pos_ce = -48:1:-36;
P_1 = -65:1:-35;
val_1 = ones(1,length(P_1));
for i = 1:length(BER_ce)
    val_1(P_1 == BER_pos_ce(i)) = BER_ce(i);
    pos_temp = find(P_1 == BER_pos_ce(i));
end

val_1(pos_temp:length(P_1)) = 0;

th = power(1 - val,200);
th_1 = power(1 - val_1,200);
%%
figure;hold on;
plot(P,th*100);
plot(P_1,th_1*100);
%%
ber_ce_pos = -57:1:-46;
ber_val_1123 = reshape(BER_CE_1123(1:36),3,12);
BER_theory_ce = 1e1*berfading(ber_ce_pos+60, 'psk', 4, 1, 100);
figure;
semilogy(ber_ce_pos,flip(ber_val_1123(1,:)) );
hold on;
semilogy(ber_ce_pos,flip(ber_val_1123(2,:)) );
semilogy(ber_ce_pos,flip(ber_val_1123(3,:)) );
semilogy(ber_ce_pos,BER_theory_ce);
%%
th1 = power(1 - flip(ber_val_1123(1,:)),20);
th2 = power(1 - flip(ber_val_1123(2,:)),20);
th_measure = power(1 - flip(ber_val_1123(3,:)),20);
figure;hold on;
plot(th1);
plot(th2);
plot(th_measure);
%%
ber_pos = -62:0.1:-51;
BER_theory = 60e1*berawgn(ber_pos+68, 'psk', 4,'nondiff');
ber_val_noCE_1123 = [5.11687175444309e-20;1.64549392387050e-16;...
    1.85979660846769e-14;4.23451965347722e-12;0.0390986485242699;...
    0.0900204240577761;0.0403400033047777;0.102857322644308;...
    0.171015895060906;0.121718666674189;0.686244516201155;1];
figure;
semilogy(ber_pos,flip(ber_val_noCE_1123));
hold on;
semilogy(ber_pos,BER_theory);

th_noCE = power(1 - flip(ber_val_noCE_1123),20);
figure; hold on;
plot(th_noCE)
plot(power(1 - BER_theory,200));

% ber value substract calculate throughput
ber_pos2 = -62:2:-51;
ber_val_noCE2 = ber_val_noCE_1123(1:2:end);
figure;
semilogy(ber_pos2,flip(ber_val_noCE2));
hold on;
semilogy(ber_pos,BER_theory);
th_noCE = power(1 - flip(ber_val_noCE2),20);
figure; hold on;
plot(ber_pos2,th_noCE);
plot(ber_pos,power(1 - BER_theory,200));
%%
figure;hold on;
plot(power(1 - BER_theory_ce,200));
plot(power(1 - BER_theory,200));
plot(th_noCE);
%%
data = convertTDMS(1,'BER_1124night.tdms');
%%
ce_file = data.Data.MeasuredData(3).Data(1:36);
ce_data = data.Data.MeasuredData(4).Data(1:36);

ber_ce_pos = -56:1:-45;
ber_ce_val = reshape(ce_data(1:36),3,12);
ber_ce_pos2 = -56:2:-45;
for i = 1:3
    ber_ce_val2(i,:) = ber_ce_val(i,2:2:end);
end

BER_theory_ce = 1e1*berfading(ber_ce_pos+60, 'psk', 4, 1, 100);
figure;
semilogy(ber_ce_pos2,flip(ber_ce_val2(1,:)) );
hold on;
semilogy(ber_ce_pos2,flip(ber_ce_val2(2,:)) );
semilogy(ber_ce_pos2,flip(ber_ce_val2(3,:)) );
semilogy(ber_ce_pos,BER_theory_ce);

th1 = power(1 - flip(ber_ce_val2(1,:)),20);
th2 = power(1 - flip(ber_ce_val2(2,:)),20);
th_measure = power(1 - flip(ber_ce_val2(3,:)),20);
figure;hold on;
% plot(ber_ce_pos2,th1);
% plot(ber_ce_pos2,th2);
plot(ber_ce_pos2,th_measure);
plot(ber_ce_pos,power(1 - BER_theory_ce,200));
%%
noce_file = data.Data.MeasuredData(3).Data(38:50);
noce_data = data.Data.MeasuredData(4).Data(38:50);
ber_noce_val = noce_data;
ber_noce_pos = -62:1:-50;
BER_theory = 80e1*berawgn(ber_noce_pos+71, 'psk', 4,'nondiff');
figure;
semilogy(ber_noce_pos,flip(noce_data));
hold on;
semilogy(ber_noce_pos,BER_theory);

th_noCE = power(1 - flip(noce_data),20);
figure;
plot(th_noCE)
%%
ber_ce_pos = -56:1:-45;
ber_ce_val = reshape(ce_data(1:36),3,12);
ber_ce_pos2 = -56:2:-45;
for i = 1:3
    ber_ce_val2(i,:) = ber_ce_val(i,1:2:end);
end
%%
ber_ce_pos = -62:0.1:-45;
BER_theory_ce = 1e1*berfading(ber_ce_pos+60, 'psk', 4, 1, 50);
% figure;
% semilogy(ber_ce_pos2,flip(ber_ce_val2(1,:)) );
% hold on; 
% semilogy(ber_ce_pos2,flip(ber_ce_val2(2,:)) );
% semilogy(ber_ce_pos2,flip(ber_ce_val2(3,:)) );
% semilogy(ber_ce_pos,BER_theory_ce);
ber_ce_val3 = flip(ber_ce_val2(3,:));
ber_ce_val3 = [1,ber_ce_val3,0];
th_measure = power(1 - ber_ce_val3,20);
th_measure_ce = pchip(-60:2:-47,th_measure,-60:1:-47);
figure;hold on;
plot((-60:1:-47)+5,th_measure_ce*100,'*','MarkerSize',10,'LineWidth',3);
plot(ber_ce_pos-0.8 + 5,100*power(1 - BER_theory_ce,200),'LineWidth',4,'Color','blue');

%%
ber_pos = -62:0.1:-50;
BER_theory = 10e1*berawgn(ber_pos+66, 'psk', 4,'nondiff');
ber_val_noCE_1123 = [5.11687175444309e-20;1.64549392387050e-16;...
    1.85979660846769e-14;4.23451965347722e-12;0.0390986485242699;...
    0.0900204240577761;0.0403400033047777;0.102857322644308;...
    0.171015895060906;0.121718666674189;0.686244516201155;1];

th_noCE = power(1 - flip(ber_val_noCE_1123),50);
% figure; hold on;
% plot(th_noCE)
% plot(power(1 - BER_theory,200));

% ber value substract calculate throughput
ber_pos2 = -62:3:-51;
ber_val_noCE2 = ber_val_noCE_1123(1:3:end);
% figure;
% semilogy(ber_pos2,flip(ber_val_noCE2));
% hold on;
% semilogy(ber_pos,BER_theory);
th_noCE2 = power(1 - flip(ber_val_noCE2),100);
th_measure_noCE = pchip(-62:3:-51,th_noCE2,-62:1:-51);
% figure; hold on;
hold on;
p2 = plot((-62:1:-51)+0.6,100*th_measure_noCE,'^','MarkerSize',10,'LineWidth',3);
% plot( -62:1:-51,th_noCE);
plot(ber_pos-0.45,100*power(1 - BER_theory,400),'LineWidth',4,'Color','red');
axis([-62 -43 0 100]);
legend([p2,p1],{'Conduted','SPAC method'});