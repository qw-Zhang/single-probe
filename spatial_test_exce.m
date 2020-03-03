%360
clear;
time = [8,9,10,12,15,18,20,30,36,40,45,72,90];
error = zeros(1,length(time));
for i = 1:length(time)
    real_phi = linspace(-pi,pi,time(i));
    [error(i), sim_sig, sim_theo] = spatial_correlation_simulation(real_phi);
%     figure;hold on;
%     plot(abs(sim_sig));
%     plot(abs(sim_theo));
end
figure;plot(error);