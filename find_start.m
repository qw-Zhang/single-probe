function location = find_start(res_corr)
len = (length(res_corr) + 1)/2;
[pks, locs] = findpeaks(abs(res_corr(len:1.5*len)));
% figure;
% plot(locs,pks);

len_pks = length(pks);
delta = 20*abs(pks(len_pks) - pks(len_pks-1));
for i = len_pks:-1:2
    temp = abs(pks(i) - pks(i - 1));
    if(temp > delta)
        break;
    end
end
location = locs(i);
end