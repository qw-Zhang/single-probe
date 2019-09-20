for i = 1:6
    for j = 1:5
        for k = 1:3
%         name_p = ['VNA_',num2str(angle_map(i)),'_90_',num2str(power_map(j)),'.dat'];
        name_p = ['VNA_',num2str(angle_map(i)),'_',num2str(power_map(j)),'_',num2str(k),'.dat'];
%         VNA_-1_90_-10_2.dat
        fid = fopen(name_p,'w+');
        end
    end
end