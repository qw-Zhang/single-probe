% real probe scenario
% different sampling function
function [real_PAS,phi_sample_res] = PAS_sample(ideal_phi,ideal_PAS,PAS_para,phi_sample,rate,sampling_method)
switch sampling_method
    case 'uniform'
        for i = 1:length(phi_sample)
            real_PAS(i) = sum(ideal_PAS(rate*(i-1)+1:rate*i))/rate;
            phi_sample_res(i) = sum(ideal_phi(rate*(i-1)+1:rate*i))/rate;
        end
    % this version just support single cluster method
    case 'narrow'
        [v, p] = max(ideal_PAS); % p is the max value position of PAS, similarly as AOA
%         AOA = PAS_para.AOA;
        AS_deg = PAS_para.AS * 180 / pi;
        % ideal_PAS times compared with 360 points.
        ideal_rate = round(length(ideal_PAS) / 360);
%         for i = 1:length(phi_sample)
%             % fixed sampling points roughly and temporary (three probes)
%             
%         end
        real_PAS(1) = ideal_PAS(p - 50);
        real_PAS(2) = ideal_PAS(p);
        real_PAS(3) = ideal_PAS(p + 50);
        
        phi_sample_res(1) = ideal_phi(p - 50);
        phi_sample_res(2) = ideal_phi(p);
        phi_sample_res(3) = ideal_phi(p + 50);
        
    otherwise
        warning('wrong sampling');
        quit;
end

real_sum_PAS = sum(real_PAS);
real_PAS = real_PAS/real_sum_PAS;

end