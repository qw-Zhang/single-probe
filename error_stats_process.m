load('D:\study\OTA\expriment_iecas\DATA\DATA0811\3D_turn\3D_turn_data_extra.mat');
% load('D:\study\OTA\expriment_iecas\DATA\DATA0811\2D_turn\turn_data_extra.mat');
% load('D:\study\OTA\expriment_iecas\DATA\DATA0810\tripod\tripod_data_extra.mat');
% load('D:\study\OTA\expriment_iecas\DATA\DATA0810\multi_tripod\multi_data_extra.mat');
% load('D:\study\OTA\expriment_iecas\DATA\DATA0824\3D_1\3D_bit_data_extra.mat')
load('D:\study\OTA\expriment_iecas\DATA\DATA0811\2D_turn\model_standard.mat');

% for c = 1:2 %amp or tau
for i = 1:6 %position
    for j = 1:4 %antenna and (correct?)
        data_error(1,i*4 + j - 4) = data_extra(1,i*4+j - 4)-model_standard(1,i);
        data_error(2,i*4 + j - 4) = (data_extra(2,i*4+j - 4)-data_extra(2,j))-model_standard(2,i)*10^6;
    end
end
% end

error_stats = zeros(2,4);
eps = 1e-10;
for j = 1:2
    for k = 1:2
        num = 0;
        for p = 1:6
            error_stats(1,2*j+k-2) = error_stats(1,2*j+k-2)+power((data_error(1,p*4+2*j+k-6)/(model_standard(1,p))),2);
%             error_stats(2,2*j+k-2) = error_stats(2,2*j+k-2)+power((data_error(2,p*4+2*j+k-6)/(model_standard(2,p)*10^6+eps)),2);
            if data_error(2,p*4+2*j+k-6)<(1/30.72)
                num = num + 1;
            end
            %                 error_stats(1,2*j+k-2) = error_stats(1,2*j+k-2)+power((data_error(1,p*4+2*j+k-6)),2);
            %                 error_stats(2,2*j+k-2) = error_stats(2,2*j+k-2)+power((data_error(2,p*4+2*j+k-6)),2);
        end
        error_stats(1,2*j+k-2) = (error_stats(1,2*j+k-2)/6);
%         error_stats(2,2*j+k-2) = (error_stats(2,2*j+k-2)/6);
        error_stats(2,2*j+k-2) = (num/6);
        %             error_stats(i,2*j+k-2) = sqrt(error_stats(i,2*j+k-2)/6);
    end
end
error_stats(1,:) = sqrt(error_stats(1,:));
