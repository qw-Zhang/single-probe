%360degree
% clear;clc;
% s = rng;

fc = 2.45e9;
c = 3e8;
lambda = c/fc;
d = linspace(0,2*lambda,41);
d_sum = linspace(0,2*lambda,41);
% probe_pos = linspace(8,36,15);  %8-36 locations of probe
probe_pos = [8,16,36,72];
error_MPAC = zeros(1,length(probe_pos));
error_SPAC = zeros(1,length(probe_pos));
phi_a = [0 30 60 120 150]*pi/180;
for j = 1:1
    for i = 3:3
        %uniform choose the location of probe
        phi_sample = linspace(-pi,pi,probe_pos(i));
%         phi_sample = linspace(-pi,pi,3); % this sample just use in spa_corr_grid_mpac_v3_1
        m = 0;v = 0;
        error_para = [m,v];
        %     [error_MPAC(i),error_SPAC(i), sim_real_sig_MPAC,sim_real_sig_SPAC, sim_sig, sim_theo,theo] = ...
        %         spa_corr_grid_simulation_v2(phi_sample,error_para,true);
        % how to implement scme CE on spatial correlation
        mpac_out = spa_corr_grid_mpac_v2_1(phi_sample,phi_a(j),d,error_para,false);
        error_MPAC(i) = mpac_out.stat;
        limits_low = [0.99 0.9 0.7 0.46 0.22 0 0 0 0 0.08 0.23];
        limits_high = [1.015 0.975 0.83 0.63 0.4 0.21 0.18 0.18 0.2 0.32 0.56];
        limits_low = [limits_low zeros(1,10)];
        limits_high = [limits_high zeros(1,10)];
        %     theo = mpac_out.theory;
        %     sim_theo = mpac_out.spatial_num;
        %     sim_real_sig_MPAC = mpac_out.spatial_circle_real;
        %     sim_ideal_sig_MPAC = mpac_out.spatial_circle;
%         figure;
%         hold on;
%         plot(d/lambda,abs(mpac_out.theory),'black');
%         plot(d/lambda,abs(mpac_out.num),'green');
% %         plot(d_sum/lambda,abs(mpac_out.spatial_num_v2),'yellow');
%         plot(d_sum/lambda,abs(mpac_out.circle_real),'red');
%         plot(d_sum/lambda,abs(mpac_out.h_sum),'magenta');
%         plot(d_sum/lambda,abs(mpac_out.h_sig),'black');
%         plot(d_sum/lambda,abs(mpac_out.h_sig_m1));
%         plot(d_sum/lambda,limits_low,'red','LineStyle','--');
%         plot(d_sum/lambda,limits_high,'red','LineStyle','-- ');
%         axis([0 2 0 1]);
%         %     plot(d/lambda,abs(sim_real_sig_SPAC),'magenta');
%         %     plot(d/lambda,abs(sim_sig),'blue');
%         xlabel('Antenna Separation in wavelength');
%         ylabel('Spatial Correlation');
%         grid on;
%         % legend('theo equation','theo num','sim circle h','sim circle sig');
% %         legend('Theory','Simulation with error position of antenna','Simulation with antenna pattern' );
%         title(num2str(phi_a(j)*180/pi));
    end
end
%
% figure;plot(probe_pos, error_SPAC);
% xlabel('number of antenna positions');
% ylabel('RMSE');