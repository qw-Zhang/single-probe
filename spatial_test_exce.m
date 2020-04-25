%360degree
clear;clc;
s = rng;
% %%
% %fix probe positions and error
% len = zeros(1,100);
% error_tmp = 1e-11*randn(1,100*length(len));
% time = linspace(1,20,20);
% error = zeros(1,length(time));
% real_phi = linspace(-pi,pi,20);
% 
% for i = 1:length(time)
%     %uniform choose the location of probe
%     error_tk = time(i)*error_tmp;
%     for j = 1:100
%         tmp = error_tk(length(real_phi)*(j-1)+1:length(real_phi)*j);
%         error_range(i,j) = max(tmp) - min(tmp);
%         error_mean(i) = mean(error_range(i,:));
%     end
%     [error(i), sim_sig, sim_theo] = spatial_correlation_simulation(real_phi,error_tk,false);
%     %     figure;hold on;
%     %     plot(abs(sim_sig));
%     %     plot(abs(sim_theo));
% end
% figure;plot(error);

%%
probe_pos = linspace(8,36,15);  %8-36 locations of probe
error_1 = zeros(1,length(probe_pos));
for i = 15:15
    %uniform choose the location of probe
    phi_sample = linspace(-pi,pi,probe_pos(i));
    error_tk = 0*randn(1,100*length(phi_sample));
    [error_1(i), sim_sig, sim_theo] = spatial_correlation_simulation_v2(phi_sample,error_tk,true);
    %     figure;hold on;
    %     plot(abs(sim_sig));
    %     plot(abs(sim_theo));
end
figure;plot(probe_pos, error_1);
% 
% %%
% probe_pos = linspace(8,36,15);  %8-36 locastions of probe
% error_1 = zeros(1,length(probe_pos));
% for i = 1:length(probe_pos)
%     %uniform choose the location of probe
%     real_phi = linspace(-pi,pi,probe_pos(i));
%     error_tk = 1e-11*randn(1,100*length(real_phi));
%     [error_1(i), sim_sig, sim_theo] = spatial_correlation_simulation(real_phi,error_tk,false);
%     %     figure;hold on;
%     %     plot(abs(sim_sig));
%     %     plot(abs(sim_theo));
% end
% figure;plot(probe_pos,error_1);
% %%
% temp = 10e-11*randn(1,10000);
% pd = fitdist(temp','Normal');
% figure;
% hist_data =  hist(temp,50);
% x = linspace(min(temp),max(temp),50);
% n = normpdf(x, pd.mu, pd.sigma);
% n = n*(max(hist_data)/max(n));
% hold on;
% plot(x,n);
% bar(x,hist_data)