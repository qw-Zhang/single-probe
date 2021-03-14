% spa_corr_grid -> this version change the way of error
% v2 -> add phase estimate
% v2.1 add simulation h(channel parameter), and simulation the spatial
% correlation without signal, just using channel itself.
% v2.3 zadoffChuSeq used
function spatial_output = spa_corr_grid_mpac_v2_3(phi_sample,phi_a,d,error_para,ant_able)
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
% t = linspace(0,1e-8,101);
% sig = exp(1i*2*pi*fc*t);
%% generate zadoffchu seq as signal
sig = zadoffChuSeq(25,101);

% figure;
% plot(abs(seq));

% seq_shift = 0.78 * seq .* exp(1i*pi/5);
% for i = 1:200
% seq_shift_awgn = awgn(seq_shift,10);
% [val(i,:),pos] = xcorr(seq_shift_awgn,seq);
% 
% ang(i) = angle(val(i,101));
% 
% end


%% PAS for theory and spatial_num value
ideal_phi = linspace(-pi,pi,3600);
% real_phi as a parameter
real_PAS = zeros(1,length(phi_sample));
rate = floor(length(ideal_phi)/length(phi_sample));

%scenario parameters
% phi_a = 0*pi/180;
scenario = 'micro';

ideal_PAS = generate_PAS1(ideal_phi,phi_sample,scenario);

% real probe scenario
for i = 1:length(phi_sample)
    real_PAS(i) = sum(ideal_PAS(rate*(i-1)+1:rate*i))/rate;
end
phi_sample = downsample(ideal_phi,rate);
real_sum_PAS = sum(real_PAS);
real_PAS = real_PAS/real_sum_PAS;

%% init matrix
h = repmat(1+1j, length(d), 2, length(ideal_phi));
[d_1, d_2, delta_d, h_sig_1, h_sig_2] = deal(zeros(length(d), length(ideal_phi)));
[d_real_1, d_real_2, delta_real_d, h_sig_real_1, h_sig_real_2,corr_h_m1,sig_err_m1]...
    = deal(zeros(length(d), length(phi_sample)));
rx_1 = repmat(1+1j,length(d),length(ideal_phi),length(sig));
rx_2 = repmat(1+1j,length(d),length(ideal_phi),length(sig));
rx_real_1 = repmat(1+1j,length(d),length(phi_sample),length(sig));
rx_real_2 = repmat(1+1j,length(d),length(phi_sample),length(sig));
[spa_sig, spa_real_sig_MPAC,spa_num,spatial,corr_h_no_delay_m2,...
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
[h_syn,h_para] = generate_H(phi_sample,scenario);
[fading_sig,h_t_tau] = generate_fading_sig_v2(h_syn,h_para,sig);

size_h = size(h_syn.h_ori);
T = size_h(4);
num_probe = size_h(3);
% connect chirp and start point
t = 0:1/1e3:1;
sig_chirp_ori = 0.5 .* chirp(t,0,1,100,'complex');
% chirp length = 1001; zeros length = 100;
sig_chirp = [zeros(1,100), sig_chirp_ori, zeros(1,100)];

% init matrix related with h matrix
[sig_err_fft_1,sig_err_fft_2,h_probe_1,h_probe_2,rx_h_sig_fft_1,rx_h_sig_fft_2] = ...
    deal(zeros(length(d),length(phi_sample),T));
[sig_err_sum_1,sig_err_sum_2] = deal(zeros(length(d),T));
% 1000 is temporary parameter.
fading_sig_part = zeros(num_probe,T,1000);
fading_sig_probe_temp = zeros(num_probe,T * 1000);
fading_sig_probe = zeros(num_probe, T * 1000 + length(sig_chirp));

for j = 1:size_h(3)
    for t_index = 1:T
        fading_sig_part(j,t_index,:) = fading_sig(j,t_index,1:1000);
    end
    fading_sig_probe_temp(j,:) = reshape(squeeze(fading_sig_part(j,:,:)).',1,T*1000);
    fading_sig_probe(j,:) = [sig_chirp, fading_sig_probe_temp(j,:)];
end

%% calculate spatial correlation value(s)
% center coordinate of two ants
%     x0 = 0.5*lambda;y0 = -0.3*lambda;
x0 = 0;y0 = 0;
%simulate radius of circle
r = 1;
beta = -pi + (pi--pi).*rand(1,length(ideal_phi)); %random select phase of signal
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
    h1 = squeeze(sum(h_syn.h1,1));
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
        
        % free space pathloss
        alpha_1 = pl_free(d_real_1(i,j),lambda);
        alpha_2 = pl_free(d_real_2(i,j),lambda);
        h_sig_real_1(i,j) = alpha*exp(1j*(beta(j) + 2*pi*fc*(d_real_1(i,j)/c) )).*sqrt(real_PAS(j));
        h_sig_real_2(i,j) = alpha*exp(1j*(beta(j) + 2*pi*fc*(d_real_2(i,j)/c) )).*sqrt(real_PAS(j));

        h_probe_1(i,j,:) = h1(j,:) .* exp(1j*k_CONST*dot(top,pos_ant_1) );
        h_probe_2(i,j,:) = h1(j,:) .* exp(1j*k_CONST*dot(top,pos_ant_2) );
        
        [ang_P_1,ang_P_2] = generate_ang_of_pattern_v2(top,pos_ant_1,pos_ant_2);
        rx_real_1(i,j,:) = P_az_amp(ang_P_1)*(h_sig_real_1(i,j)*sig);
        rx_real_2(i,j,:) = P_az_amp(ang_P_2)*(h_sig_real_2(i,j)*sig);
        
        rx_fad_sig_probe_1(j,:) = fading_sig_probe(j,:) .* exp(1j*k_CONST*dot(top,pos_ant_1));
        rx_fad_sig_probe_2(j,:) = fading_sig_probe(j,:) .* exp(1j*k_CONST*dot(top,pos_ant_2));
        
        % no process of chirp
        % that means bias one or two samples would influence the res of
        % xcorr
        xcorr_rx_1(i,j,:) = xcorr_process(rx_fad_sig_probe_1(j,:), sig, sig_chirp_ori, T);
        xcorr_rx_2(i,j,:) = xcorr_process(rx_fad_sig_probe_2(j,:), sig, sig_chirp_ori, T);
        
        xcorr_m1_temp(i,j) = mean(xcorr_rx_1(i,j,:) .* conj(xcorr_rx_2(i,j,: )));
        
        for t_index = 1:T
            rx_h_sig_1 = fading_sig(j,t_index,1:1000).*exp(1j*k_CONST*dot(top,pos_ant_1));
            rx_h_sig_2 = fading_sig(j,t_index,1:1000).*exp(1j*k_CONST*dot(top,pos_ant_2));
            rx_h_sig_fft_1(i,j,t_index) = fft_process( rx_h_sig_1 );
            rx_h_sig_fft_2(i,j,t_index) = fft_process( rx_h_sig_2 );
            %                 tx = awgn(sig,10);
            %                 tx_z = awgn(sig_z,10);
            
            if error_para(2) ~= 0
                % connect a sin (single tone)
                lf = 50;
                t_sin = 0:0.001:0.05;
                temp_sin = exp(1i*pi*lf*t_sin);
                % this version suppose that the delay estimation is
                % precise.
                rx_sig_err_temp_1 = [temp_sin,squeeze(fading_sig(j,t_index,1:1000))'];
                rx_sig_err_temp_1 = rx_sig_err_temp_1 .* exp(1j*k_CONST*dot(top,pos_ant_1));
                rx_sig_err_temp_2 = [temp_sin,squeeze(fading_sig(j,t_index,1:1000))'];
                rx_sig_err_temp_2 = rx_sig_err_temp_2 .* exp(1j*k_CONST*dot(top,pos_ant_2));
                
                sample_t = 400;
                
                [ang_est,error] = deal(zeros(1,sample_t));
                
                for k = 1:sample_t
                    % error of time alignment on probe
                    error(k) = error_para(1) + error_para(2)*randn(1,1);
                    rx_sig_err_2_sample = rx_sig_err_temp_2 .* exp(1i*2*pi*fc*error(k));
                    % angle estimate just use the single tone signal
                    fft_temp = fft(rx_sig_err_2_sample(1:length(temp_sin)));
                    [v_max, pos_max] = max(abs(fft_temp));
                    ang_est(k) = angle(fft_temp(pos_max));
                end
                
                [h_y,h_x] = hist(ang_est);   %choose a most probably piont by histogram
                [vh_y,ph_y] = max(h_y);
                %decrease origin phase and using new phase
                % rx_real_2(i,j,:) = rx_real_sample_2*exp(1i*-ang_est(k))*exp(1i*h_x(ph_y));
                phase_est = h_x(ph_y);
                % phase_est = mean(ang_est);
                len_rx_sig_err_temp_1 = length(rx_sig_err_temp_1);
                rx_sig_err_1 = rx_sig_err_temp_1(length(temp_sin)+1:len_rx_sig_err_temp_1);
                
                len_rx_sig_err_temp_2 = length(rx_sig_err_temp_2);
                rx_sig_est_2 = rx_sig_err_2_sample(length(temp_sin)+1:len_rx_sig_err_temp_2)...
                               *exp(1i*-ang_est(k))*exp(1i*phase_est);
                
                sig_err_fft_1(i,j,t_index) = fft_process( rx_sig_err_1 );
                sig_err_fft_2(i,j,t_index) = fft_process( rx_sig_est_2 );
            end
            
        end
        
        sig_err_m1(i,j) = mean(sig_err_fft_1(i,j,:) .* conj(sig_err_fft_2(i,j,:)));
        
        fad_sig_m1_temp(j) = mean(rx_h_sig_fft_1(i,j,:) .* conj(rx_h_sig_fft_2(i,j,:)));
        
        corr_h_m1_temp = corrcoef(h_probe_1(i,j,:),h_probe_2(i,j,:));
        corr_h_m1(i,j) = corr_h_m1_temp(1,2);
    end
    
    corr_fad_sig_m1(i) = sum(fad_sig_m1_temp);
    
    corr_fad_sig_err_m1(i) = sum(sig_err_m1(i,:));
    
    sig_err_sum_1(i,:) = squeeze(sum(sig_err_fft_1(i,:,:),2));
    sig_err_sum_2(i,:) = squeeze(sum(sig_err_fft_2(i,:,:),2));
    corr_fad_sig_err_m2(i) = mean(sig_err_sum_1(i,:) .* conj(sig_err_sum_2(i,:)));
    
    sig_sum_1 = squeeze(sum(rx_h_sig_fft_1,2));
    sig_sum_2 = squeeze(sum(rx_h_sig_fft_2,2));
    fad_sig_m2_temp = corrcoef(sig_sum_1(i,:),sig_sum_2(i,:));
    corr_fad_sig_m2(i) = fad_sig_m2_temp(1,2);
    
    xcorr_m1(i) = sum(xcorr_m1_temp(i,:));
    
    xcorr_rx_sum_1 = squeeze(sum(xcorr_rx_1(i,:,:),2));
    xcorr_rx_sum_2 = squeeze(sum(xcorr_rx_2(i,:,:),2));
    xcorr_m2(i) = mean(xcorr_rx_sum_1 .* conj(xcorr_rx_sum_2));
    
    corr_para_bakup(i) = sum(sig_sum_1(:,i) .* conj(sig_sum_2(:,i)));
    
    % no time delay in h
    h_1(i,:) = squeeze(sum(h_probe_1(i,:,:),2))';
    h_2(i,:) = squeeze(sum(h_probe_2(i,:,:),2))';
    corr_h_no_delay_m2(i) = mean(h_1(i,:) .* conj(h_2(i,:)));
    
    %randomly choose one point in signal sequence. "100" is chosen randomly
    num = randi([1,length(sig)],1);
    spa_real_sig_MPAC(i) = sum(rx_real_1(i,:,num).*conj(rx_real_2(i,:,num)));
    
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
corr_fad_sig_err_m1 = corr_fad_sig_err_m1 ./ corr_fad_sig_err_m1(1);
corr_fad_sig_err_m2 = corr_fad_sig_err_m2 ./ corr_fad_sig_err_m2(1);
xcorr_m1 = xcorr_m1 ./ xcorr_m1(1);
xcorr_m2 = xcorr_m2 ./ xcorr_m2(1);

%% plot figure
figure;
hold on;
plot(d/lambda,abs(spatial),'black');
plot(d/lambda,abs(spa_num),'green');
%         plot(d_sum/lambda,abs(mpac_out.spatial_num_v2),'yellow');
plot(d/lambda,abs(spa_real_sig_MPAC),'red');
% plot(d/lambda,abs(spatial_output.h_sum),'magenta');
plot(d/lambda,abs(corr_fad_sig_m2),'black');
plot(d/lambda,abs(corr_fad_sig_m1),'magenta');
plot(d/lambda,abs(corr_fad_sig_err_m1),'blue');
plot(d/lambda,abs(corr_fad_sig_err_m2),'cyan');
plot(d/lambda,abs(xcorr_m1));
plot(d/lambda,abs(xcorr_m2));
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
