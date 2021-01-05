%spa_corr_grid -> this version change the way of error
%v2 -> add phase estimate
function spatial_output = spa_corr_grid_mpac_v2(phi_sample,phi_a,d,error_para,ant_able)
    fc = 2.45e9;
    c = 3e8;
    lambda = c/fc;
%     new_d = linspace(0,lambda,100);
%     d = linspace(0,lambda,length(new_d));
    new_d = linspace(0,d(length(d))/2,length(d));
%     d = linspace(0,new_d(length(new_d)),length(new_d));
    ideal_phi = linspace(-pi,pi,3600);
    
    % real_phi as a parameter
    % real_phi = linspace(-pi,pi,18);
    real_PAS = zeros(1,length(phi_sample));
    
    rate = floor(length(ideal_phi)/length(phi_sample));
    % rate = length(ideal_phi)/length(real_phi);
    
    %scenario parameters
%     phi_a = 0*pi/180;
    % choose one scenario from "test" "micro" "macro"
    scenario = 'micro';
    ideal_PAS = generate_PAS1(ideal_phi,scenario);
    
    % real probe scenario
    for i = 1:length(phi_sample)
        real_PAS(i) = sum(ideal_PAS(rate*(i-1)+1:rate*i))/rate;
        phi_sample(i) = sum(ideal_phi(rate*(i-1)+1:rate*i))/rate;
    end
    real_sum_PAS = sum(real_PAS);
    real_PAS = real_PAS/real_sum_PAS;
    
    %Rx one antenna pattern
    %generate the pattern gain of the antenna.
    %parameter:antenna,frequency,al,ez.
    %G(dBi)=10lgGi G(dBd)=10lgG dBi=dBd+2.15
    % [D,al,ez] = pattern(ant,fc,0:1:360,-90:1:90);
    
    if ant_able
        P_az = generate_pattern(fc);
    else
        P_az = zeros(1,3600);
    end
    P_az_amp = 10.^(P_az./10);
    
    spatial = zeros(1,length(new_d));
%     beta = -pi + (pi--pi).*rand(1,length(ideal_phi)); %random select phase of signal
    beta = zeros(1,length(ideal_phi));
    alpha = 1;
    
    %using numberical method calculate the spatial correlation
    spatial_num = zeros(length(new_d),1);
    
    %simulate radius of circle
    r = 1;
    
    t = linspace(0,1e-9,100);
    t = t(1:50);
    % doppler frequency fd = fc*v*cos(theta)/c
    v = 10; %m/s
    theta = 30*pi/180;
    fd = fc*v*cos(theta)/c;
    phi_sig = -pi + 2*pi*randn(1,length(t));
%     sig = exp(1i*phi_sig);
    sig = exp(1i*2*pi*(fc)*t);
    
    %init matrix
    h = repmat(1+1j, length(new_d), 2, length(ideal_phi));
    [d_1, d_2, delta_d, h_sig_1, h_sig_2] = deal(zeros(length(new_d), length(ideal_phi)));
    [d_real_1, d_real_2, delta_real_d, h_sig_real_1, h_sig_real_2] = deal(zeros(length(new_d), length(phi_sample)));
    rx_1 = repmat(1+1j,length(new_d),length(ideal_phi),length(sig));
    rx_2 = repmat(1+1j,length(new_d),length(ideal_phi),length(sig));
    rx_real_1 = repmat(1+1j,length(new_d),length(phi_sample),length(sig));
    rx_real_2 = repmat(1+1j,length(new_d),length(phi_sample),length(sig));
    [spatial_circle, spatial_circle_sig, spatial_circle_real_sig_MPAC,...
        spatial_circle_real] = deal(zeros(1,length(new_d)));
    h_mpac = zeros(2,length(new_d));
    [rx_1_sum,rx_2_sum] = deal(zeros(length(new_d),length(sig)));
    % snr = 20;
%     x0 = 0.5*lambda;y0 = -0.3*lambda;   %center coordinate of two ants
    x0 = 0;y0 = 0;
    
    for i = 1:length(new_d)
        %simulate spatial correlation using two antennas with circle
        pos_ant_1 = [new_d(i)*cos(phi_a),new_d(i)*sin(phi_a)]; %(x,y)
        pos_ant_2 = [-new_d(i)*cos(phi_a),-new_d(i)*sin(phi_a)];
        
        ang_P = zeros(2,length(ideal_phi));
        for j = 1:length(ideal_phi)
            phi = ideal_phi(j) + pi/2;
%             top = [r*cos(phi+pi/2),r*sin(phi+pi/2)];
            top = [r*cos(phi),r*sin(phi)];
            d_1(i,j) = sqrt((pos_ant_1(1) - top(1))^2 + (pos_ant_1(2) - top(2))^2);
            d_2(i,j) = sqrt((pos_ant_2(1) - top(1))^2 + (pos_ant_2(2) - top(2))^2);
            [ang_P(1,j),ang_P(2,j)] = generate_ang_of_pattern_v2(top,pos_ant_1,pos_ant_2);
            delta_d(i,j) = d_2(i,j) - d_1(i,j);
            h_sig_1(i,j) = alpha*exp(1j*(2*pi*(dot(pos_ant_1,top))/lambda)).*sqrt(ideal_PAS(j));
            h_sig_2(i,j) = alpha*exp(1j*(2*pi*(dot(pos_ant_2,top))/lambda)).*sqrt(ideal_PAS(j));
%             [ang_P_1,ang_P_2] = generate_ang_of_pattern(d_1(i,j),d_2(i,j),new_d(i),r,(pi/2 - phi),(phi + pi/2));
            rx_1(i,j,:) = P_az_amp(ang_P(1,j))*h_sig_1(i,j)*sig;
            rx_2(i,j,:) = P_az_amp(ang_P(2,j))*h_sig_2(i,j)*sig;
        end
        
        rx_1_sum(i,:) = reshape(sum(rx_1(i,:,:),2),1,length(sig));
        rx_2_sum(i,:) = reshape(sum(rx_2(i,:,:),2),1,length(sig));
        % Corr(:,:,i) = abs(corrcoef(h_new_1(i,:),h_new_2(i,:)));
        % spatial_circle_sig(i) = sum(r1(i,:).*conj(r2(i,:)))/sum(r1(1,:).*conj(r2(1,:)));
        
        num = randi([1,length(sig)],1);
        %randomly choose one point in signal sequence. "100" is chosen randomly
        spatial_circle_sig(i) = sum(rx_1(i,:,num).*conj(rx_2(i,:,num)));
        %     spatial_circle(i) = sum(h_sig_1(i,:).*conj(h_sig_2(i,:)));
        
        %simulate spatial correlation using two antenna with circle(MPAC),
        %real phi and PAS
        for j = 1:length(phi_sample)
            phi_real = phi_sample(j) + pi/2;
            %%postioner
            error_top = 0.03*randn*lambda;
%             error_top = 0;
            top = [error_top + r*cos(phi_real),error_top + r*sin(phi_real)];   %this errro is used to first positioner
            pos_ant_1 = [x0 + new_d(i)*cos(phi_a),y0 + new_d(i)*sin(phi_a)]; %(x,y)
            pos_ant_2 = [x0 - new_d(i)*cos(phi_a),y0 - new_d(i)*sin(phi_a)];
            d_real_1(i,j) = sqrt((pos_ant_1(1) - top(1))^2 + (pos_ant_1(2) - top(2))^2);
            d_real_2(i,j) = sqrt((pos_ant_2(1) - top(1))^2 + (pos_ant_2(2) - top(2))^2);
            delta_real_d(i,j) = d_real_2(i,j) - d_real_1(i,j);
            
            %error of time alignment on probe
            error = error_para(1) + error_para(2)*randn(1,1);
            % ideal equal pathloss
            % beta angle shows the separate the measurement like SPAC.
            % In MPAC, beta is the different distance of multiple probes to
            % DUT antennas, it rely on the phi_sample (j).
            h_sig_real_1(i,j) = alpha*exp(1j*(beta(j) + 2*pi*fc*(d_real_1(i,j)/c) )).*sqrt(real_PAS(j));
            h_sig_real_2(i,j) = alpha*exp(1j*(beta(j) + 2*pi*fc*(d_real_2(i,j)/c) )).*sqrt(real_PAS(j));
            [ang_P_1(j),ang_P_2(j)] = generate_ang_of_pattern_v2(top,pos_ant_1,pos_ant_2);
            rx_real_1(i,j,:) = P_az_amp(ang_P_1(j))*(h_sig_real_1(i,j)*sig);
            rx_real_2(i,j,:) = P_az_amp(ang_P_2(j))*(h_sig_real_2(i,j)*sig);
            
            % using single DUT antenna method()
            
        end
        
        rx_real_1_sum(i,:) = reshape(sum(rx_real_1(i,:,:),2),1,length(sig));
        rx_real_2_sum(i,:) = reshape(sum(rx_real_2(i,:,:),2),1,length(sig));
        
        %randomly choose one point in signal sequence. "100" is chosen randomly
        
        %%num = randi([1,length(sig)],1);
        spatial_circle_real_sig_MPAC(i) = sum(rx_real_1(i,:,num).*conj(rx_real_2(i,:,num)));
        
        
        beta1 = 0;
        %using numberical method calculate the spatial correlation
        for j = 1:length(ideal_phi)
            h(i,1,j) = alpha*exp(1j*(beta1)).*sqrt(ideal_PAS(j));
            h(i,2,j) = alpha*exp(1j*(beta1 + 2*pi*(2*new_d(i)/lambda) * cos((pi/2-ideal_phi(j)) + phi_a) ) ).*sqrt(ideal_PAS(j));
        end
        h1_cal = reshape(h(i,1,:),1,length(ideal_PAS));
        h2_cal = reshape(conj(h(i,2,:)),1,length(ideal_PAS));
        spatial_num(i,:) = sum(h1_cal.*h2_cal);
        
        %     Corr(i,:,:) = abs(corrcoef(squeeze(h(i,1,:)),squeeze(h(i,2,:))));
        
        %using traditional method(ideally equation) calculate
        tau = zeros(1,length(ideal_phi));
        for k = 1:length(ideal_phi)
            tau(k) = d(i)*sin(ideal_phi(k)-phi_a)/c;
            %                 spatial(2,i) = spatial(2,i) + exp(-1i*2*pi*fc * (tau(k) + error_para(k)) ).*PAS(k);
            spatial(i) = spatial(i) + exp(-1i*2*pi*fc * tau(k) ).*(ideal_PAS(k));
        end
        ss(i) = mean(rx_real_1_sum(i,:).*conj(rx_real_2_sum(i,:)));
    end
    spatial = spatial./spatial(1);
    spatial_circle_sig = spatial_circle_sig./spatial_circle_sig(1);
    spatial_circle_real_sig_MPAC = spatial_circle_real_sig_MPAC./spatial_circle_real_sig_MPAC(1);
    % Correaltion = squeeze(Corr(:,1,2));
    
    %calculate the statistical error
    stat_MPAC = sqrt(sum(power(abs(spatial_circle_real_sig_MPAC)-abs(spatial_circle_sig),2))/length(spatial_circle_sig));
    ss = ss./ss(1);
%     output structure
spatial_output.stat = stat_MPAC;
spatial_output.spatial_circle_real = spatial_circle_real_sig_MPAC;
spatial_output.spatial_circle = spatial_circle_sig;
spatial_output.spatial_num = spatial_num;
spatial_output.theory = spatial;
spatial_output.spatial_sum = ss;

% different strength of signal to control probe power
