com_tx = Txdata(1,1:50485) + 1i*Txdata(2,1:50485);

% temp = data110dbm2(1,1:300000) + 1i*data110dbm2(2,1:300000);
for i = 1:6
    temp_ch = num2str(i);
    var_name(i,:) = ['data',temp_ch,'10dbm2'];
end
%%
% figure;hold on;
for i = 1:6
    temp = eval(var_name(i,:));
    data(i,:) = temp(1,1:300000) + 1i*temp(2,1:300000);
    [a,b] = xcorr(data(i,:), com_tx);
    [v,p] = max(abs(a));
    res(i,2) = v;
    res(i,1) = b(p);
%     plot(b,abs(a));
end

%%
% 频域计算GCC-PHAT
Nfft = 1024*300;
% figure;hold on;
for i = 1:1
%     temp = eval(var_name(i,:));
%     data = temp(1,1:300000) + 1i*temp(2,1:300000);
    Gss = fft(data(i,:),Nfft).*conj(fft(com_tx,Nfft));
    Gss = Gss./abs(Gss);
    % Gss = exp(1i*angle(Gss)); % 这种方式也可以
    xcorrGcc = fftshift(ifft(Gss));
    xaxis = linspace(-Nfft/2,Nfft/2,Nfft);
    [v,p] = max(abs(xcorrGcc));
    res(i,4) = v;
    res(i,3) = xaxis(p);
    %     plot(xaxis,abs(xcorrGcc));
%     tt = Txdata(1,1:300000) + 1i*Txdata(2,1:300000);
%     [aa,bb]=gccphat(data(i,:),tt);
end
%%
for j = 1:2        
    tv_res = res(1,j * 2);
    tp_res = res(1,j * 2 - 1);
    for i = 1:6
        res(i,j * 2) = res(i,j * 2)/tv_res;
        res(i,j * 2 - 1) = res(i,j * 2 - 1) - tp_res;
    end
end
%%
%cc(sum)
data_sum = sum(data,1);
    [a,b] = xcorr(data_sum, com_tx);
    figure;plot(b,abs(a));
    
%%
% gcc(sum)
data_sum = sum(data,1);
Gss = fft(data_sum,Nfft).*conj(fft(com_tx,Nfft));
Gss = Gss./abs(Gss);
% Gss = exp(1i*angle(Gss)); % 这种方式也可以
xcorrGcc = fftshift(ifft(Gss));
xaxis = linspace(-Nfft/2,Nfft/2,Nfft);
figure;plot(xaxis,abs(xcorrGcc));
%%
temp_sum_res = gcc;
for i = 1:6
    res_sum(i,2) = temp_sum_res(i).Position(2)/temp_sum_res(1).Position(2);
end