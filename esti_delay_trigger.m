clear;
% test = load('..\DATA0115\test.dat');

file_path ='D:\study\OTA\expriment_iecas\DATA\DATA0119';
data_list=dir(strcat(file_path,'\*.dat'));
len_dir = length(data_list);

loc_start = zeros(len_dir,2);
for j = 1:1:len_dir
    data_name = data_list(j).name;
    test=load(data_name);
    
    n = 1;
    t = linspace(0,n,n*100);
    tx = sin(2*pi*1*t)+1i*sin(2*pi*1*t);
    
    rx(1,:) = test(1,:)+1i*test(2,:);
    rx(2,:) = test(3,:)+1i*test(4,:);
    len = length(rx(1,:));
    sig = zeros(2,len);
    for i = 1:1:2
        repair = IQrepair(rx(i,:),tx);
        if repair
            sig(i,:) = -1.*real(rx(i,:)) + 1i.*imag(rx(i,:));
        else
            sig(i,:) = rx(i,:);
        end
        
        [a,b] = xcorr(sig(i,:),tx);
        % figure;
        % plot(b,a);
        loc_start(j,i) = find_start(a) - 10;
    end
end

%%
