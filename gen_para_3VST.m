%%
ang = [1 -13 146 -30 -11 -1];
ang = deg2rad(ang);
com = cos(ang) + 1i*sin(ang);
ang_std = [0 120 240];
ang_std = deg2rad(ang_std);
com_std = cos(ang_std) + 1i*sin(ang_std);

for i = 1:6
%     temp_ch = num2str(i);
%     name(i,:) = ['case',temp_ch,'.dat'];
%     temp_tx = load(name(i,:));
%     tx(i,:) = temp_tx(1,:)+1i*temp_tx(2,:);
%     for j = 1:3
% %         tx_de(i*3+j-3,:) = tx(i,:).*cos(ang(i)-ang_std(j));
%             tx_angle(i*3+j-3) = cos(ang(i)-ang_std(j));
%     end
%         
%     
rr(:,i) = gen_para_3(ang(i),i);
end

%%
result_tx = zeros(3,length(tx));
result_txangle = zeros(3,1);
cchirp = real(tx(1,100001:110000));
for i = 1:3
    for j = 1:6
        result_tx(i,:) = result_tx(i,:)+tx(j,:).*rr(i,j);
%             result_txangle(i) = result_txangle(i) + tx_angle(3*j+i-3);
    end
    
    result_ttx(2*i-1,:) = real(result_tx(i,:));
    result_ttx(2*i-1,100001:110000) = cchirp;
    result_ttx(2*i,:) = imag(result_tx(i,:));
    result_ttx(2*i,100001:110000) = cchirp;
end


    
%%
fid = fopen('Txdata_3VSt.dat','w+');
for i =1:6
    for j = 1:length(result_ttx(i,:))
        fprintf(fid,'%f\t',result_ttx(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);

