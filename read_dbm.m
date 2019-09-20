data = textread('VST0_0_10.txt','%s','delimiter','\n');
data_1 = data(275);
data_2 = cell2mat(deblank(data_1));
data_3 = str2num(data_2);
dbm = data_3(2)
v = sqrt((50*power(10,(dbm/10)))/1000)