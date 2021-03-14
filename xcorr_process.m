function [res,xcorr_res] = xcorr_process(rx_sig,ori_sig,preamble,T)
%% find start point
[val_pre, pos_pre] = xcorr(rx_sig,preamble);
[v_max_pre,p_max_pre] = max(abs(val_pre));
start_pos = pos_pre(p_max_pre);
start_pos = start_pos + 1101;
% judge a point
%% cross-correlation process
batch_size = T;
each_size = 1000;

res = zeros(1,T);
xcorr_res = zeros(T,each_size * 2 - 1);

for i = 1:batch_size
    data_temp = rx_sig( (start_pos + (i - 1) * each_size + 1 ) : (start_pos + i * each_size) );
    val = xcorr(data_temp,ori_sig');
    [v_max,p_max] = max(abs(val));
    res(i) = val(p_max);
    xcorr_res(i,:) = val;
end