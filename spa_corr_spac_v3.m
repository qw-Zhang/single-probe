%test with h1 h2 without the euqalation of saptial correlation
%using three method generate spatial correlation:
%   1. traditional formula method
%   2. numberical calculation method
%   3. simulation two antennas with circle method
%   4. add signal transfer based method_3

%this version change the way of error
function [stat, spatial_circle_real_sig,spatial_circle_sig, spatial_num] = ...
    spa_corr_spac_v3(phi_sample,error_para,ant_able)
fc = 2.535e9;   %carrier frequency
c = 3e8;        %speed of microwave
lambda = c/fc;  %wavelength
d = linspace(0,lambda,100); %array of distance in two antennas
new_d = linspace(0,lambda/2,100);   %array of distance in one antenna with origin
ideal_phi = linspace(-pi,pi,3600);  %ideal division of 360 degrees

% real_phi as a parameter
% real_phi = linspace(-pi,pi,18);
real_PAS = zeros(1,length(phi_sample)); %division number of actual angle with phi_sample

rate = floor(length(ideal_phi)/length(phi_sample)); %
% rate = length(ideal_phi)/length(real_phi);

%scenario parameters
phi_a = 0*pi/180;   %angle of init arrangment (the center of the line between of two antennas)
scenario = 'micro'; %choose scenario(below three scenarios)
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
        error('wrong scenario.');
end

PAS = generate_PAS(ideal_phi,pow,AS,step_mov);  %generate ideal PAS(length of ideal_PAS)


% real probe scenario
% extract using real PASs from ideal PAS and sum up them
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

% using special antenna's pattern 
if ant_able
    P_az = generate_pattern(fc);
else
    P_az = zeros(1,360);
end
P_az_amp = 10.^(P_az./10);

spatial = zeros(2,length(d));
beta = -pi + (pi--pi).*rand(1,length(ideal_phi)); %random select phase of signal
alpha = 1;  % fading factor of signal

%using numberical method calculate the spatial correlation
spatial_num = zeros(length(d),1);

%simulate radius of circle
r = 1;

t = linspace(0,1e-9,100);
% doppler frequency fd = fc*v*cos(theta)/c
v = 10; %m/s
theta = 30*pi/180;  % direction of movement
fd = fc*v*cos(theta)/c; %generate doppler frequency
sig = exp(1i*2*pi*(fc + fd)*t); %transmitted siganl

%init matrix
h = repmat(1+1j, length(d), 2, length(ideal_phi));
[d_1, d_2, delta_d, h_sig_1, h_sig_2] = deal(zeros(length(d), length(ideal_phi)));
[d_real_1, d_real_2, delta_real_d, h_sig_real_1, h_sig_real_2] = deal(zeros(length(d), length(phi_sample)));
rx_1 = repmat(1+1j,length(d),length(ideal_phi),length(sig));
rx_2 = repmat(1+1j,length(d),length(ideal_phi),length(sig));
rx_real_1 = repmat(1+1j,length(d),length(phi_sample),length(sig));
rx_real_2 = repmat(1+1j,length(d),length(phi_sample),length(sig));
sample_t = 1000;    %sample times at one angle
error = zeros(1,sample_t);
ang_error = zeros(1,sample_t);
error_mean = zeros(length(d),length(phi_sample));
% rx_real_sample_1 = repmat(1+1j,length(d),length(phi_sample),length(sig));
% rx_real_sample_2 = repmat(1+1j,sample_t,length(sig));
rx_real_sample_2 = repmat(1+1j,1,length(sig));
[spatial_circle, spatial_circle_sig, spatial_circle_real_sig, spatial_circle_real] = deal(zeros(1,length(d)));
ang_est = zeros(1,sample_t);
ang_ori = zeros(length(d),length(phi_sample));
ang_debug = zeros(length(phi_sample),2);

for i = 1:length(d)
    %simulate spatial correlation using two antennas with circle
    %this method using the geometric relationship
    
    % ideal arrival angle of center
    for j = 1:length(ideal_phi)
        phi = ideal_phi(j) - phi_a; %arrival angle this tap
        %distance from transmitter to antenna 1
        d_1(i,j) = sqrt( r^2 + new_d(i)^2 - 2*new_d(i)*r*cos( pi/2 - phi ) );
        %distance from transmitter to antenna 2
        d_2(i,j) = sqrt( r^2 + new_d(i)^2 - 2*new_d(i)*r*cos( phi + pi/2 ) );
        delta_d(i,j) = d_2(i,j) - d_1(i,j);
        %channel factor for antenna 1 and 2 with fading factor, phase and PAS
        h_sig_1(i,j) = alpha*exp(1j*(beta(j))).*sqrt(PAS(j));
        h_sig_2(i,j) = alpha*exp(1j*(beta(j) + 2*pi*delta_d(i,j)/lambda)).*sqrt(PAS(j));
        %generate pattern factor corresponding the situation and
        %convolution the signal
        [ang_P_1,ang_P_2] = generate_ang_of_pattern(d_1(i,j),d_2(i,j),new_d(i),r,(pi/2 - phi),(phi + pi/2));
        rx_1(i,j,:) = P_az_amp((ang_P_1))*conv(h_sig_1(i,j),sig);
        rx_2(i,j,:) = P_az_amp((ang_P_2))*conv(h_sig_2(i,j),sig);
    end
    
    % this method is from 977, but it doesn't work
    % Corr(:,:,i) = abs(corrcoef(h_new_1(i,:),h_new_2(i,:)));
    % spatial_circle_sig(i) = sum(r1(i,:).*conj(r2(i,:)))/sum(r1(1,:).*conj(r2(1,:)));
    
    num = randi([1,length(d)],1);
    %randomly choose one point in signal sequence. "100" is chosen randomly
    spatial_circle_sig(i) = sum(rx_1(i,:,num).*conj(rx_2(i,:,num)));
    %generate spatial correlation ideally
    spatial_circle(i) = sum(h_sig_1(i,:).*conj(h_sig_2(i,:)));
    
    %simulate spatial correlation using two antennas with circle,
    %real phi and PAS
    
    % this block is same as above block, just using real locations
    for j = 1:length(phi_sample)
        phi_real = (phi_sample(j) - phi_a);
        %         phi_real_2 = pi - phi_real_1;
        d_real_1(i,j) = sqrt( r^2 + new_d(i)^2 - 2*new_d(i)*r*cos( pi/2 - phi_real ) );
        d_real_2(i,j) = sqrt( r^2 + new_d(i)^2 - 2*new_d(i)*r*cos( phi_real + pi/2 ) );
        delta_real_d(i,j) = d_real_2(i,j) - d_real_1(i,j);
        h_sig_real_1(i,j) = alpha*exp(1j*(beta(j))).*sqrt(real_PAS(j));
        h_sig_real_2(i,j) = alpha*exp(1j*(beta(j) + 2*pi*delta_real_d(i,j)/lambda)).*sqrt(real_PAS(j));
        [ang_P_1,ang_P_2] = generate_ang_of_pattern(d_real_1(i,j),d_real_2(i,j),new_d(i),r,(pi/2 - phi_real),(phi_real + pi/2));
%         ang_debug(j,:) = [ang_P_1, ang_P_2];
        rx_real_1(i,j,:) = P_az_amp(ang_P_1)*(h_sig_real_1(i,j)*sig);
        rx_real_2(i,j,:) = P_az_amp(ang_P_2)*(h_sig_real_2(i,j)*sig);
%         sig_ori_fft = fft(rx_real_2(i,j,:),512);
%         [v_ori_max, pos_ori_max] = max(abs(sig_ori_fft));
%         ang_ori(i,j) = angle(sig_ori_fft(pos_ori_max));

        %this loop use selecting best sample for eliminate the influence of
        %synchornation with two recievers
        for k = 1:sample_t
            error(k) = error_para(1) + error_para(2)*randn(1,1); %genenrate error
            rx_real_sample_2 = rx_real_2(i,j,:)*exp(1i*2*pi*fc*error(k));   %signal sample
            fft_temp = fft(rx_real_sample_2,512);   %generate phase by fft
            [v_max, pos_max] = max(abs(fft_temp));  %choose the max piont
            ang_est(k) = angle(fft_temp(pos_max));  %get the every phase 
        end
        [h_y,h_x] = hist(ang_est,50);   %choose a most probably piont by histogram
        [vh_y,ph_y] = max(h_y);
        %decrease origin phase and using new phase
        rx_real_2(i,j,:) = rx_real_sample_2*exp(1i*-ang_est(k))*exp(1i*h_x(ph_y));
    end
    
    %randomly choose one point in signal sequence. "100" is chosen randomly
    %     for k = 1:length(real_phi)
    %         num = randi([1,length(sig)],1);
    %         spatial_circle_real_sig(i) = spatial_circle_real_sig(i) + (rx_real_1(i,k,num).*conj(rx_real_2(i,k,num)));
    %     end
    num = randi([1,length(sig)],1);
    %generate spatial correlation
    spatial_circle_real_sig(i) = sum(rx_real_1(i,:,num).*conj(rx_real_2(i,:,num)));
    
    %     spatial_circle_real_sig_w2(i) =
    spatial_circle_real(i) = sum(h_sig_real_1(i,:).*conj(h_sig_real_2(i,:)));
    
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
    
%     %using traditional method(ideally equation) calculate
%     tau = zeros(1,length(ideal_phi));
%     for k = 1:length(ideal_phi)
%         tau(k) = d(i)*sin(ideal_phi(k)-phi_a)/c;
%         spatial(2,i) = spatial(2,i) + exp(-1i*2*pi*fc * (tau(k) + error_para(k)) ).*PAS(k);
%     end
end
spatial_circle_sig = spatial_circle_sig./spatial_circle_sig(1);
spatial_circle_real_sig = spatial_circle_real_sig./spatial_circle_real_sig(1);
% Correaltion = squeeze(Corr(:,1,2));

%calculate the statistical error
stat = sqrt(sum(power(abs(spatial_circle_real_sig)-abs(spatial_circle_sig),2))/length(spatial_circle_sig));

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
% legend('theo num','sim circle real sig','sim circle sig');