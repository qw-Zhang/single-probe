% generate PAS data by function:PAS_fun
function [ideal_PAS, real_PAS] = generate_PAS1(ideal_phi,phi_sample,scenario)
    switch scenario
        case 'test'
            AOA = -22.5*pi/180;AS = 35*pi/180;
            pow = 0;

        case 'single'
            AOA = 0*pi/180;AS = 35*pi/180;
            pow = 0;

        case 'micro'
            % % SCME urban micro-cell
            AOA = [0.7,-13.2,146.1,-30.5,-11.4,-1.1]*pi/180;AS = 35*pi/180;
            pow = [-3,-4.3,-5.7,-7.3,-9,-11.4];

        case 'macro'
            % %SCME urban macro-cell
            AOA = [65.7,45.6,143.2,32.5,-91.1,-19.2]*pi/180;AS = 35*pi/180;
            pow = [-3,-5.2,-4.7,-8.2,-12.1,-15.5];

        otherwise
            warning('wrong scenario.');
            quit;
    end

    amp = 10.^(pow/10);
    A = amp./max(amp);
    PAS_fun = @(x) (1/(2*AS))*sum( A.*exp(-abs(x-AOA)/AS) );

    PAS = zeros(1,length(ideal_phi));
    for i = 1:length(ideal_phi)
        PAS(i) = PAS_fun(ideal_phi(i));
    end
    sumPAS = sum(PAS);
    ideal_PAS = PAS./sumPAS;
    PAS_para.AOA = AOA;
    PAS_para.AS = AS;

end

% %%
% % generate PAS data by function:PAS_fun
% function PAS_output = generate_PAS_inline(phi_sample,scenario,cluster_index)
% switch scenario
%     case 'test'
%         AOA = -22.5*pi/180;AS = 35*pi/180;
%         pow = 0;
%         
%     case 'micro'
%         % % SCME urban micro-cell
%         AOA = [0.7,-13.2,146.1,-30.5,-11.4,-1.1] * pi/180;
%         %         AOA = [0.7,-13.2,146.1,-30.5,-21.4,30.1] * pi/180;
%         AS = 35*pi/180;
%         pow_table = [ [-3 -5.2 -7]; [-4.3 -6.5 -8.3]; [-5.7 -7.9 -9.7]; [-7.3 -9.5 -11.3]; [-9 -11.2 -13]; [-11.4 -13.6 -15.4] ];
%         pow = mean(pow_table,2);
%         
%     case 'macro'
%         % %SCME urban macro-cell
%         AOA = [65.7,45.6,143.2,32.5,-91.1,-19.2]*pi/180;
%         AS = 35*pi/180;
%         pow = [-3,-5.2,-4.7,-8.2,-12.1,-15.5];
%         
%     otherwise
%         warning('wrong scenario.');
%         quit;
% end
% 
% if strcmp(scenario, 'test')
%     cluster_index = 1;
% end
% 
% % PAS_fun = @(x) (1/(2*AS))* exp(-abs(x - AOA(cluster_index) )/AS);
% PAS_fun = @(x) (1/(2*AS))* exp(-abs(x - 0 )/AS);
% 
% ideal_phi = linspace(-pi,pi,3600);
% ideal_PAS_temp = zeros(size(ideal_phi));
% 
% for ideal_index = 1:length(ideal_phi)
%     ideal_PAS_temp(ideal_index) = PAS_fun( ideal_phi(ideal_index) );
% end
% 
% move_step = round(AOA(cluster_index) * 180 / pi * 10);
% ideal_PAS = circshift(ideal_PAS_temp,move_step);
% ideal_PAS = ideal_PAS ./ max(ideal_PAS);
% rate = floor(length(ideal_phi)/length(phi_sample));
% % real probe scenario
% for i = 1:length(phi_sample)
%     PAS_output.PAS(i) = sum(ideal_PAS(rate*(i-1)+1:rate*i))/rate;
%     phi_sample(i) = sum(ideal_phi(rate*(i-1)+1:rate*i))/rate;
% end
% real_PAS_temp = downsample(ideal_PAS,rate);
% PAS_output.phi_sample = downsample(ideal_phi,rate);
% real_sum_PAS = sum(real_PAS_temp);
% % real_PAS = real_PAS_temp./real_sum_PAS;
% real_PAS = real_PAS_temp;
% 
% % just for test
% PAS_output.ideal_PAS = ideal_PAS;
% PAS_output.pow_init = pow(1);
% % PAS_output.PAS = real_PAS;
% PAS_output.pow = pow(cluster_index);
% end