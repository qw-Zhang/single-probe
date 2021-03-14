% spa_corr_grid -> this version change the way of error
% v2 -> add phase estimate
% v2.1 add simulation h(channel parameter), and simulation the spatial
% correlation without signal, just using channel itself.
% v2.4 delete fading H and sig, just save (ori) synthesized PAS
function spatial_output = spa_corr_grid_mpac_v2_4(phi_sample,phi_a,d,error_para,ant_able)
%% some basic parameters
fc = 2.45e9;
c = 3e8;
speed_of_light=2.99792458e8;
wavelength=speed_of_light/fc;
k_CONST = 2*pi/wavelength;
lambda = c/fc;

new_d = linspace(0,d(length(d))/2,length(d));

% Fading info are involved in H matrix.
% sig is a standard single tone signal.
t = linspace(0,1e-8,101);
sig = exp(1i*2*pi*fc*t);

%% PAS for theory and spatial_num value
ideal_phi = linspace(-pi,pi,3600);
% real_phi as a parameter
real_PAS = zeros(1,length(phi_sample));
rate = floor(length(ideal_phi)/length(phi_sample));

%scenario parameters
% phi_a = 0*pi/180;
scenario = 'micro';

ideal_PAS = generate_PAS1(ideal_phi,phi_sample,scenario);
PP = load('ideal_PP_micro.mat');
ideal_PAS = PP.ideal_PP_micro;
% real probe scenario
for i = 1:length(phi_sample)
    real_PAS(i) = sum(ideal_PAS(rate*(i-1)+1:rate*i))/rate;
end
phi_sample = downsample(ideal_phi,rate);
real_sum_PAS = sum(real_PAS);
real_PAS = real_PAS/real_sum_PAS;
% macro
% real_PAS = [0.0701    0.0277    0.0414    0.0431    0.1138    0.3280    0.2040    0.1720];
% micro
real_PAS = [0.2150    0.1633    0.2401    0.4448    0.6384    0.3471    0.2322    0.2993];
%% init matrix
h = repmat(1+1j, length(d), 2, length(ideal_phi));
[d_1, d_2, delta_d, h_sig_1, h_sig_2] = deal(zeros(length(d), length(ideal_phi)));
[d_real_1, d_real_2, delta_real_d, h_sig_real_1, h_sig_real_2,corr_h_m1,sig_err_m1]...
    = deal(zeros(length(d), length(phi_sample)));
rx_1 = repmat(1+1j,length(d),length(ideal_phi),length(sig));
rx_2 = repmat(1+1j,length(d),length(ideal_phi),length(sig));
rx_real_1 = repmat(1+1j,length(d),length(phi_sample),length(sig));
rx_real_2 = repmat(1+1j,length(d),length(phi_sample),length(sig));
[spa_sig, spa_real_sig_MPAC, spa_real_sig_MPAC_m2,spa_num,spatial,corr_h_no_delay_m2,...
    corr_fad_sig_m1,corr_fad_sig_m2,corr_para_bakup] = deal(zeros(1,length(d)));

%% Rx one antenna pattern
%generate the pattern gain of the antenna.
%parameter:antenna,frequency,al,ez.
%G(dBi)=10lgGi G(dBd)=10lgG dBi=dBd+2.15
if ant_able
    P_az = generate_pattern(fc);
else
    P_az = zeros(1,3600);
end
P_az_amp = 10.^(P_az./20);

%% generate H fading matrix and fading signal

%% calculate spatial correlation value(s)
% center coordinate of two ants
%     x0 = 0.5*lambda;y0 = -0.3*lambda;
x0 = 0;y0 = 0;
%simulate radius of circle
r = 1;
beta = -pi + (pi--pi).*rand(1,length(ideal_phi)); %random select phase of signal
beta = pi / 4 * ones(1,length(ideal_phi));
alpha = 1;
for i = 1:length(new_d)
    %simulate spatial correlation using two antennas with circle
    pos_ant_1 = [new_d(i)*cos(phi_a),new_d(i)*sin(phi_a)]; %(x,y)
    pos_ant_2 = [-new_d(i)*cos(phi_a),-new_d(i)*sin(phi_a)];
    for j = 1:length(ideal_phi)
        phi = ideal_phi(j) - phi_a;
        top = [r*cos(phi+pi/2),r*sin(phi+pi/2)];
        d_1(i,j) = sqrt((pos_ant_1(1) - top(1))^2 + (pos_ant_1(2) - top(2))^2);
        d_2(i,j) = sqrt((pos_ant_2(1) - top(1))^2 + (pos_ant_2(2) - top(2))^2);
        delta_d(i,j) = d_2(i,j) - d_1(i,j);
        h_sig_1(i,j) = alpha*exp(1j*(beta(j) + 2*pi*dot(pos_ant_1,top)/lambda)).*sqrt(ideal_PAS(j));
        h_sig_2(i,j) = alpha*exp(1j*(beta(j) + 2*pi*dot(pos_ant_2,top)/lambda)).*sqrt(ideal_PAS(j));
        %             [ang_P_1,ang_P_2] = generate_ang_of_pattern(d_1(i,j),d_2(i,j),new_d(i),r,(pi/2 - phi),(phi + pi/2));
        [ang_P_1,ang_P_2] = generate_ang_of_pattern_v2(pos_ant_1,pos_ant_2,top);
        rx_1(i,j,:) = P_az_amp((ang_P_1))*conv(h_sig_1(i,j),sig);
        rx_2(i,j,:) = P_az_amp((ang_P_2))*conv(h_sig_2(i,j),sig);
    end
    
    num = randi([1,length(d)],1);
    %randomly choose one point in signal sequence. "100" is chosen randomly
    spa_sig(i) = sum(rx_1(i,:,num).*conj(rx_2(i,:,num)));
    %     spatial_circle(i) = sum(h_sig_1(i,:).*conj(h_sig_2(i,:)));
    
    %simulate spatial correlation using two antenna with circle(MPAC),
    %real phi and PAS
%     h1 = squeeze(sum(h_syn.h1,1));
    for j = 1:length(phi_sample)
        phi_real = phi_sample(j) + pi/2;
        %%postioner error
        %         error_top = 0.3*randn*lambda;
        error_top = 0;
        top = [error_top + r*cos(phi_real),error_top + r*sin(phi_real)];   %this errro is used to first positioner
        pos_ant_1 = [x0 + new_d(i)*cos(phi_a),y0 + new_d(i)*sin(phi_a)]; %(x,y)
        pos_ant_2 = [x0 - new_d(i)*cos(phi_a),y0 - new_d(i)*sin(phi_a)];
        d_real_1(i,j) = sqrt((pos_ant_1(1) - top(1))^2 + (pos_ant_1(2) - top(2))^2);
        d_real_2(i,j) = sqrt((pos_ant_2(1) - top(1))^2 + (pos_ant_2(2) - top(2))^2);
        delta_real_d(i,j) = d_real_2(i,j) - d_real_1(i,j);
        
        h_sig_real_1(i,j) = alpha*exp(1j*(beta(j) + 2*pi*fc*(d_real_1(i,j)/c) )).*sqrt(real_PAS(j));
        h_sig_real_2(i,j) = alpha*exp(1j*(beta(j) + 2*pi*fc*(d_real_2(i,j)/c) )).*sqrt(real_PAS(j));

%         h_probe_1(i,j,:) = h1(j,:) .* exp(1j*k_CONST*dot(top,pos_ant_1) );
%         h_probe_2(i,j,:) = h1(j,:) .* exp(1j*k_CONST*dot(top,pos_ant_2) );
        
        [ang_P_1,ang_P_2] = generate_ang_of_pattern_v2(top,pos_ant_1,pos_ant_2);
        rx_real_1(i,j,:) = P_az_amp(ang_P_1)*(h_sig_real_1(i,j)*sig);
        rx_real_2(i,j,:) = P_az_amp(ang_P_2)*(h_sig_real_2(i,j)*sig);
        
    end
    
    %randomly choose one point in signal sequence. "100" is chosen randomly
%     num = randi([1,length(sig)],1);
    num = 10;
    spa_real_sig_MPAC(i) = sum(rx_real_1(i,:,num).*conj(rx_real_2(i,:,num)));
    rx_real_sum_1(i) = sum(rx_real_1(i,:,num),2);
    rx_real_sum_2(i) = sum(rx_real_2(i,:,num),2);
    spa_real_sig_MPAC_m2(i) = sum(rx_real_sum_1(i).*conj(rx_real_sum_2(i)));
    
    %using numberical method calculate the spatial correlation
    for j = 1:length(ideal_phi)
        h(i,1,j) = alpha*exp(1j*(beta(j))).*sqrt(ideal_PAS(j));
        h(i,2,j) = alpha*exp(1j*(beta(j) + 2*pi*(2*new_d(i)/lambda) * cos((pi/2-ideal_phi(j)) + phi_a) ) ).*sqrt(ideal_PAS(j));
    end
    h1_cal = reshape(h(i,1,:),1,length(ideal_PAS));
    h2_cal = reshape(conj(h(i,2,:)),1,length(ideal_PAS));
    spa_num(i,:) = sum(h1_cal.*h2_cal);
    
    %     Corr(i,:,:) = abs(corrcoef(squeeze(h(i,1,:)),squeeze(h(i,2,:))));
    
    %using traditional method(ideally equation) calculate
    tau = zeros(1,length(ideal_phi));
    for k = 1:length(ideal_phi)
        tau(k) = d(i)*sin(ideal_phi(k)-phi_a) / c;
        spatial(i) = spatial(i) + exp(-1i*2*pi*fc * tau(k) ) .* ideal_PAS(k);
    end
end

%% normalize correlation
corr_h_no_delay_m2 = corr_h_no_delay_m2 ./ corr_h_no_delay_m2(1);
corr_fad_sig_m1 = corr_fad_sig_m1 ./ corr_fad_sig_m1(1);
spatial = spatial./spatial(1);
spa_sig = spa_sig./spa_sig(1);
spa_real_sig_MPAC = spa_real_sig_MPAC./spa_real_sig_MPAC(1);
spa_real_sig_MPAC_m2 = spa_real_sig_MPAC_m2 ./ spa_real_sig_MPAC_m2(1);
%% plot figure
figure;
hold on;
plot(d/lambda,abs(spatial),'black');
plot(d/lambda,abs(spa_num),'green');
%         plot(d_sum/lambda,abs(mpac_out.spatial_num_v2),'yellow');
plot(d/lambda,abs(spa_real_sig_MPAC),'red');
plot(d/lambda,abs(spa_real_sig_MPAC_m2),'blue');
axis([0 2 0 1]);
%     plot(d/lambda,abs(sim_real_sig_SPAC),'magenta');
%     plot(d/lambda,abs(sim_sig),'blue');

if isequal(scenario,'micro')
    limits_low = [0.99 0.9 0.7 0.46 0.22 0 0 0 0 0.08 0.23];
    limits_high = [1.015 0.975 0.83 0.63 0.4 0.21 0.18 0.18 0.2 0.32 0.56];
    remain = length(d) - length(limits_low);
    limits_low = [limits_low zeros(1,remain)];
    limits_high = [limits_high zeros(1,remain)];
    plot(d/lambda,limits_low,'red','LineStyle','--');
    plot(d/lambda,limits_high,'red','LineStyle','--');
end

xlabel('Antenna Separation in wavelength');
ylabel('Spatial Correlation');
grid on;

%% return structure
%calculate the statistical error
stat_MPAC = sqrt(sum(power(abs(spa_real_sig_MPAC)-abs(spa_sig),2))/length(spa_sig));

spatial_output.stat = stat_MPAC;
spatial_output.h_sum = corr_h_no_delay_m2;
spatial_output.h_sig = corr_fad_sig_m2;
spatial_output.h_sig_m1 = corr_fad_sig_m1;
spatial_output.circle_real = spa_real_sig_MPAC;
spatial_output.circle = spa_sig;
spatial_output.num = spa_num;
spatial_output.theory = spatial;

end
