% h = zeros(S,CLUSTER,K,T); resolution:10ns
% sig(t) = linspace(0,1e-8,101); resolution:10ns
%
% 1. generate delay domain(tau) vector, length decided by delta_t (like scme_test.m)
% 2. convolution the h(tau) with sig
% 3. connection the convoloed signal by time domain (reshape)
function [fading_sig,h_t_tau] = generate_fading_sig(h,h_para,k_index,sig)
    h_ori = h.h_ori;
    h1 = h.h1;
    
    delay = h_para.delay;
    delta_T = h_para.delta_T;
    size_h = size(h_ori);
    T = size_h(4);
%     K = size_h(3);
    CLUSTER = size_h(2);
    
    % resolution of sig: 1e-8
    sig_T = length(sig) * 1e-8;
    vec_tau_len = round(delta_T / sig_T);
    vec_tau = zeros(1,vec_tau_len);
    
    for clu_index = 1:CLUSTER
        for t_index = 1:T
            vec_tau(t_index,delay(clu_index)) = h1(clu_index,k_index,t_index);
        end
    end

    
    for t_index = 1:T
        sig_temp = conv(squeeze(vec_tau(t_index, :)),sig);
        fading_sig(t_index,:) = sig_temp;
    end    
    
    h_t_tau = vec_tau;
end