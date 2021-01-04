% real probe scenario
% different sampling function
function real_PAS = PAS_sample(ideal_phi,phi_sample,rate,sampling_method)
switch sampling_method
    case 'uniform'
        for i = 1:length(phi_sample)
            real_PAS(i) = sum(ideal_PAS(rate*(i-1)+1:rate*i))/rate;
            phi_sample(i) = sum(ideal_phi(rate*(i-1)+1:rate*i))/rate;
        end
        
    case 'narrow'
        for i = 1:length(phi_sample)
            
        
    otherwise
        warning('wrong sampling');
        quit;
end

real_sum_PAS = sum(real_PAS);
real_PAS = real_PAS/real_sum_PAS;

end