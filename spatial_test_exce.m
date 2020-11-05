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
fc = 2.535e9;
c = 3e8;
lambda = c/fc;
d = linspace(0,lambda,100);
% probe_pos = linspace(8,36,15);  %8-36 locations of probe
probe_pos = [8,16,36,72];
error_1 = zeros(1,length(probe_pos));

for i = 4:4
    %uniform choose the location of probe
    phi_sample = linspace(-pi,pi,probe_pos(i));
    error_tk = 0*randn(1,100*length(phi_sample));
    m = 0;v = 0;
    error_para = [m,v];
    [error_1(i), sim_real_sig, sim_sig, sim_theo] = spa_corr_grid_simulation_v2(phi_sample,error_para,true);
figure;
hold on;
% plot(d/lambda,abs(spatial(2,:)),'black');
plot(d/lambda,abs(sim_theo),'green');
plot(d/lambda,abs(sim_real_sig),'red');
plot(d/lambda,abs(sim_sig),'blue');
xlabel('Antenna Separation in wavelength');
ylabel('Spatial Correlation');
grid on;
% legend('theo equation','theo num','sim circle h','sim circle sig');
legend('theo num','sim circle real sig','sim circle sig' );
end

figure;plot(probe_pos, error_1);
xlabel('number of antenna positions');
ylabel('RMSE');