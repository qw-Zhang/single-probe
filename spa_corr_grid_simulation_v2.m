%spa_corr_grid -> this version change the way of error
%v2 -> add phase estimate
function [stat_MPAC,stat_SPAC, spatial_circle_real_sig_MPAC,spatial_circle_real_sig_SPAC,spatial_circle_sig, spatial_num,spatial] = ...
        spa_corr_grid_simulation_v2(phi_sample,error_para,ant_able)
    fc = 2.535e9;
    c = 3e8;
    lambda = c/fc;
    d = linspace(0,lambda,100);
    new_d = linspace(0,lambda/2,100);
    ideal_phi = linspace(-pi,pi,3600);
    
    % real_phi as a parameter
    % real_phi = linspace(-pi,pi,18);
    real_PAS = zeros(1,length(phi_sample));
    
    rate = floor(length(ideal_phi)/length(phi_sample));
    % rate = length(ideal_phi)/length(real_phi);
    
    %scenario parameters
    phi_a = 0*pi/180;
    scenario = 'micro';
    switch scenario
        case 'test'
            AOA = -22.5*pi/180;AS = 35*pi/180;
            step_mov = floor(-22.5);
            pow = 0;
            
        case 'micro'
            % % SCME urban micro-cell
            AOA = [0.7,-13.2,146.1,-30.5,-11.4,-1.1]*pi/180;AS = 35*pi/180;
            step_mov = floor([0.7,-13.2,146.1,-30.5,-11.4,-1.1]);
            pow = [-3,-4.3,-5.7,-7.3,-9,-11.4];
            
        case 'macro'
            % %SCME urban macro-cell
            AOA = [65.7,45.6,143.2,32.5,-91.1,-19.2]*pi/180;AS = 35*pi/180;
            step_mov = floor([65.7,45.6,143.2,32.5,-91.1,-19.2]);
            pow = [-3,-5.2,-4.7,-8.2,-12.1,-15.5];
            
        otherwise
            warning('wrong scenario.');
            quit;
    end
    
    PAS = generate_PAS(ideal_phi,pow,AS,step_mov);
    
    
    % real probe scenario
    for i = 1:length(phi_sample)
        real_PAS(i) = sum(PAS(rate*(i-1)+1:rate*i))/rate;
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
        P_az = zeros(1,360);
    end
    P_az_amp = 10.^(P_az./10);
    
    spatial = zeros(2,length(d));
    beta = -pi + (pi--pi).*rand(1,length(ideal_phi)); %random select phase of signal
    alpha = 1;
    
    %using numberical method calculate the spatial correlation
    spatial_num = zeros(length(d),1);
    
    %simulate radius of circle
    r = 1;
    
    t = linspace(0,1e-7,100);
    % doppler frequency fd = fc*v*cos(theta)/c
    v = 10; %m/s
    theta = 30*pi/180;
    fd = fc*v*cos(theta)/c;
    sig = exp(1i*2*pi*(fc + fd)*t);
    
    %init matrix
    h = repmat(1+1j, length(d), 2, length(ideal_phi));
    [d_1, d_2, delta_d, h_sig_1, h_sig_2] = deal(zeros(length(d), length(ideal_phi)));
    [d_real_1, d_real_2, delta_real_d, h_sig_real_1, h_sig_real_2] = deal(zeros(length(d), length(phi_sample)));
    rx_1 = repmat(1+1j,length(d),length(ideal_phi),length(sig));
    rx_2 = repmat(1+1j,length(d),length(ideal_phi),length(sig));
    rx_real_1 = repmat(1+1j,length(d),length(phi_sample),length(sig));
    rx_real_2 = repmat(1+1j,length(d),length(phi_sample),length(sig));
    [spatial_circle, spatial_circle_sig, spatial_circle_real_sig_MPAC,...
        spatial_circle_real_sig_SPAC, spatial_circle_real] = deal(zeros(1,length(d)));
    ang_debug = zeros(length(phi_sample),2);
    sample_t = 2000;
    [error,ang_est] = deal(zeros(1,sample_t));
    % snr = 20;
%     x0 = 0.5*lambda;y0 = -0.3*lambda;   %center coordinate of two ants
    x0 = 0;y0 = 0;
    
    for i = 1:length(new_d)
        %simulate spatial correlation using two antennas with circle
        pos_ant_1 = [new_d(i)*cos(phi_a),new_d(i)*sin(phi_a)]; %(x,y)
        pos_ant_2 = [-new_d(i)*cos(phi_a),-1*new_d(i)*sin(phi_a)];
        for j = 1:length(ideal_phi)
            phi = ideal_phi(j) - phi_a;
            top = [r*cos(phi+pi/2),r*sin(phi+pi/2)];
            d_1(i,j) = sqrt((pos_ant_1(1) - top(1))^2 + (pos_ant_1(2) - top(2))^2);
            d_2(i,j) = sqrt((pos_ant_2(1) - top(1))^2 + (pos_ant_2(2) - top(2))^2);
            delta_d(i,j) = d_2(i,j) - d_1(i,j);
            h_sig_1(i,j) = alpha*exp(1j*(beta(j))).*sqrt(PAS(j));
            h_sig_2(i,j) = alpha*exp(1j*(beta(j) + 2*pi*delta_d(i,j)/lambda)).*sqrt(PAS(j));
            [ang_P_1,ang_P_2] = generate_ang_of_pattern(d_1(i,j),d_2(i,j),new_d(i),r,(pi/2 - phi),(phi + pi/2));
            rx_1(i,j,:) = P_az_amp((ang_P_1))*conv(h_sig_1(i,j),sig);
            rx_2(i,j,:) = P_az_amp((ang_P_2))*conv(h_sig_2(i,j),sig);
        end
        
        % Corr(:,:,i) = abs(corrcoef(h_new_1(i,:),h_new_2(i,:)));
        % spatial_circle_sig(i) = sum(r1(i,:).*conj(r2(i,:)))/sum(r1(1,:).*conj(r2(1,:)));
        
        num = randi([1,length(d)],1);
        %randomly choose one point in signal sequence. "100" is chosen randomly
        spatial_circle_sig(i) = sum(rx_1(i,:,num).*conj(rx_2(i,:,num)));
        %     spatial_circle(i) = sum(h_sig_1(i,:).*conj(h_sig_2(i,:)));
        
        %simulate spatial correlation using two antennas with circle(SPAC),
        %real phi and PAS
        
        for j = 1:length(phi_sample)
            phi_real = (phi_sample(j) - phi_a);
            %%top = [2.5*lambda + r*cos(phi_real+pi/2),-1.2*lambda +
            %%r*sin(phi_real+pi/2)];     %this error is used to second
            %%postioner
            error_top = 0.3*randn*lambda;
%             error_top = 0;
            top = [error_top + r*cos(phi_real+pi/2),error_top + r*sin(phi_real+pi/2)];   %this errro is used to first positioner
            pos_ant_1 = [x0 + new_d(i)*cos(phi_a),y0 + new_d(i)*sin(phi_a)]; %(x,y)
            pos_ant_2 = [x0 - new_d(i)*cos(phi_a),y0 - new_d(i)*sin(phi_a)];
            d_real_1(i,j) = sqrt((pos_ant_1(1) - top(1))^2 + (pos_ant_1(2) - top(2))^2);
            d_real_2(i,j) = sqrt((pos_ant_2(1) - top(1))^2 + (pos_ant_2(2) - top(2))^2);
            delta_real_d(i,j) = d_real_2(i,j) - d_real_1(i,j);
            h_sig_real_1(i,j) = alpha*exp(1j*(beta(j))).*sqrt(real_PAS(j));
            h_sig_real_2(i,j) = alpha*exp(1j*(beta(j) + 2*pi*delta_real_d(i,j)/lambda)).*sqrt(real_PAS(j));
            [ang_P_1,ang_P_2] = generate_ang_of_pattern_v2(top,pos_ant_1,pos_ant_2);
            %         ang_debug(j,:) = [ang_P_1, ang_P_2];
            rx_real_1(i,j,:) = P_az_amp(ang_P_1)*(h_sig_real_1(i,j)*sig);
            rx_real_2(i,j,:) = P_az_amp(ang_P_2)*(h_sig_real_2(i,j)*sig);
            %         rx_real_1(i,j,:) = awgn(rx_real_1(i,j,:),snr);
            %         rx_real_2(i,j,:) = awgn(rx_real_2(i,j,:),snr);
            if error_para(2) ~= 0
                for k = 1:sample_t
                    error(k) = error_para(1) + error_para(2)*randn(1,1);
                    rx_real_sample_2 = rx_real_2(i,j,:)*exp(1i*2*pi*fc*error(k));
                    fft_temp = fft(rx_real_sample_2,512);
                    [v_max, pos_max] = max(abs(fft_temp));
                    ang_est(k) = angle(fft_temp(pos_max));
                end
                [h_y,h_x] = hist(ang_est,50);
                [vh_y,ph_y] = max(h_y);
                rx_real_2(i,j,:) = rx_real_sample_2*exp(1i*-ang_est(k))*exp(1i*h_x(ph_y));
            end
        end
        
        %randomly choose one point in signal sequence. "100" is chosen randomly
        %     for k = 1:length(real_phi)
        %         num = randi([1,length(sig)],1);
        %         spatial_circle_real_sig(i) = spatial_circle_real_sig(i) + (rx_real_1(i,k,num).*conj(rx_real_2(i,k,num)));
        %     end
        num = randi([1,length(sig)],1);
        spatial_circle_real_sig_SPAC(i) = sum(rx_real_1(i,:,num).*conj(rx_real_2(i,:,num)));
        
        %     spatial_circle_real_sig_w2(i) =
        spatial_circle_real(i) = sum(h_sig_real_1(i,:).*conj(h_sig_real_2(i,:)));
        
        %simulate spatial correlation using two antenna with circle(MPAC),
        %real phi and PAS
        for j = 1:length(phi_sample)
            phi_real = (phi_sample(j) - phi_a);
            %%postioner
            error_top = 0.3*randn*lambda;
%             error_top = 0;
            top = [error_top + r*cos(phi_real+pi/2),error_top + r*sin(phi_real+pi/2)];   %this errro is used to first positioner
            pos_ant_1 = [x0 + new_d(i)*cos(phi_a),y0 + new_d(i)*sin(phi_a)]; %(x,y)
            pos_ant_2 = [x0 - new_d(i)*cos(phi_a),y0 - new_d(i)*sin(phi_a)];
            d_real_1(i,j) = sqrt((pos_ant_1(1) - top(1))^2 + (pos_ant_1(2) - top(2))^2);
            d_real_2(i,j) = sqrt((pos_ant_2(1) - top(1))^2 + (pos_ant_2(2) - top(2))^2);
            delta_real_d(i,j) = d_real_2(i,j) - d_real_1(i,j);
            
            %error of time alignment on probe
            error = error_para(1) + error_para(2)*randn(1,1);
            
            h_sig_real_1(i,j) = alpha*exp(1j*(beta(j) + 2*pi*fc*error)).*sqrt(real_PAS(j));
            h_sig_real_2(i,j) = alpha*exp(1j*(beta(j) + 2*pi*delta_real_d(i,j)/lambda + 2*pi*fc*error)).*sqrt(real_PAS(j));
            [ang_P_1,ang_P_2] = generate_ang_of_pattern_v2(top,pos_ant_1,pos_ant_2);
            rx_real_1(i,j,:) = P_az_amp(ang_P_1)*(h_sig_real_1(i,j)*sig);
            rx_real_2(i,j,:) = P_az_amp(ang_P_2)*(h_sig_real_2(i,j)*sig);
        end
        
        %randomly choose one point in signal sequence. "100" is chosen randomly
        num = randi([1,length(sig)],1);
        spatial_circle_real_sig_MPAC(i) = sum(rx_real_1(i,:,num).*conj(rx_real_2(i,:,num)));
        
        
        
        beta1 = 0;
        %using numberical method calculate the spatial correlation
        for j = 1:length(ideal_phi)
            h(i,1,j) = alpha*exp(1j*(beta1)).*sqrt(PAS(j));
            h(i,2,j) = alpha*exp(1j*(beta1 + 2*pi*(2*new_d(i)/lambda) * cos((pi/2-ideal_phi(j)) + phi_a) ) ).*sqrt(PAS(j));
        end
        h1_cal = reshape(h(i,1,:),1,length(PAS));
        h2_cal = reshape(conj(h(i,2,:)),1,length(PAS));
        spatial_num(i,:) = sum(h1_cal.*h2_cal);
        
        %     Corr(i,:,:) = abs(corrcoef(squeeze(h(i,1,:)),squeeze(h(i,2,:))));
        
            %using traditional method(ideally equation) calculate
            tau = zeros(1,length(ideal_phi));
            for k = 1:length(ideal_phi)
                tau(k) = d(i)*sin(ideal_phi(k)-phi_a)/c;
%                 spatial(2,i) = spatial(2,i) + exp(-1i*2*pi*fc * (tau(k) + error_para(k)) ).*PAS(k);
                spatial(2,i) = spatial(2,i) + exp(-1i*2*pi*fc * (tau(k)) ).*PAS(k);
            end
    end
    spatial_circle_sig = spatial_circle_sig./spatial_circle_sig(1);
    spatial_circle_real_sig_SPAC = spatial_circle_real_sig_SPAC./spatial_circle_real_sig_SPAC(1);
    spatial_circle_real_sig_MPAC = spatial_circle_real_sig_MPAC./spatial_circle_real_sig_MPAC(1);
    % Correaltion = squeeze(Corr(:,1,2));
    
    %calculate the statistical error
    stat_SPAC = sqrt(sum(power(abs(spatial_circle_real_sig_SPAC)-abs(spatial_circle_sig),2))/length(spatial_circle_sig));
    stat_MPAC = sqrt(sum(power(abs(spatial_circle_real_sig_MPAC)-abs(spatial_circle_sig),2))/length(spatial_circle_sig));
    
    % %plot spatial correlation
    % figure;
    % %plot(d/lambda,abs(spatial(1,:)),'b');
    % hold on;
    % % plot(d/lambda,abs(spatial(2,:)),'black');
    % plot(d/lambda,abs(spatial_num),'green');
    % plot(d/lambda,abs(spatial_circle_real_sig),'red');
    % plot(d/lambda,abs(spatial_circle_sig),'blue');
    % % plot(d/lambda,Correaltion,'p');
    % % axis([0 1 0 1]);
    % xlabel('Antenna Separation in wavelength');
    % ylabel('Spatial Correlation');
    % grid on;
    % % legend('theo equation','theo num','sim circle h','sim circle sig');
    % legend('theo num','sim circle real sig' );