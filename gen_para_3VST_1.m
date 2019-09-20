

%%
ang = [1 -13 146 -30 -11 -1];
ang = deg2rad(ang);
for i = 1:6
    temp = num2str(i);
    name(i,:) = ['case',temp,'.dat'];
%     load(name);
    
    var_name(i,:) = ['case',temp];
    temp = eval(var_name(i,:));
    com(i,:) = temp(1,:).*cos(ang(i)) + 1i*temp(2,:).*sin(ang(i));
    
end

%%
for i = 1:3
    res_com(i,:) = com(2*i-1,:)+com(2*i,:);
    result(2*i-1,:) = real(res_com(i,:));
    result(2*i,:) = imag(res_com(i,:));
end

%%
fid = fopen('txdata_3.dat','w+');
for i =1:6
    for j = 1:length(result(i,:))
        fprintf(fid,'%f\t',result(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);
%%
% amp = [1 0.5346 0.7471 0.3718 0.2504 0.1435];
com = amp.*cos(ang) + 1i*amp.*sin(ang);
for i = 1:3
    res_com(i) = com(2*i-1)+com(2*i);
end

amp_para = abs(res_com);
ang_para = rad2deg(angle(res_com));

Txdata_com = Txdata(1,:) + 1i*Txdata(2,:);
