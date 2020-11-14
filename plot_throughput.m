load('BER_data.mat');
P = -64:2:-36;
val_ce = zeros(1,length(P));
val = zeros(1,length(P));
BER_pos = -60:2:-28;

for i = 1:length(BER_val)
    val(P == BER_pos(i)) = 1 - BER_val(i);
end

BER_pos_ce = -48:1:-40;
P_1 = -65:1:-35;
val_1 = zeros(1,length(P_1));
for i = 1:length(BER_ce)
    val_1(P_1 == BER_pos_ce(i)) = (1 - BER_ce(i));
    pos_temp = find(P_1 == BER_pos_ce(i));
end
val_1(pos_temp:length(P_1)) = 1;

P_new = -65:0.1:-35;
val_sp = pchip(P,val,P_new);
val_ce_sp = pchip(P_1,val_1,P_new);

figure;hold on;
plot(P,val*100,'Marker','o','LineWidth',2);
plot(P_new,val_sp*100,'LineWidth',3);
plot(P_1,val_1*100,'Marker','*','LineWidth',2);
plot(P_new,val_ce_sp*100,'LineWidth',3);