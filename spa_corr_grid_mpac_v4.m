% spac_corr_grid_simulation_mpac_v3 -> two antennas change to single ref ant
% v4, add CE influnce on signal. If it works, this h parameter would
% implemented on VST CE.
function spatial_output = spa_corr_grid_simulation_mpac_v4(phi_sample,phi_a,d,error_para,ant_able)
fc = 2.45e9;
c = 3e8;
lambda = c/fc;
% new_d = linspace(0,2*lambda,100);
% d = new_d;
new_d = d;
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
    P_az = zeros(1,length(ideal_phi));
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

t = linspace(0,1e-8,500);
sig = exp(1i*2*pi*(fc)*t);
% doppler frequency fd = fc*v*cos(theta)/c
% the influnce of doppler frequency has implemented on h(channel parameter)
%     v = 10; %m/s
% v = 30*1000/(60*60);    % 30km/hour
% theta = 30*pi/180;
% fd = fc*v*cos(theta)/c;
% phi_sig = -pi + 2*pi*randn(1,length(t));
%     sig = exp(1i*phi_sig);
%% init matrix
h = repmat(1+1j, length(new_d), 2, length(ideal_phi));
[d_1, d_2, delta_d, h_sig_1, h_sig_2] = deal(zeros(length(new_d), length(ideal_phi)));
[d_real_1, d_real_2, delta_real_d, h_sig_real_1, h_sig_real_2] = deal(zeros(length(new_d), length(phi_sample)));
rx_1 = repmat(1+1j,length(new_d),length(ideal_phi),length(sig));
rx_2 = repmat(1+1j,length(new_d),length(ideal_phi),length(sig));
rx_real_1 = repmat(1+1j,length(new_d),length(phi_sample),length(sig));
rx_real_2 = repmat(1+1j,length(new_d),length(phi_sample),length(sig));
%     [spatial_circle, spatial_circle_sig, spatial_circle_real_sig_MPAC,...
%         spatial_circle_real] = deal(zeros(1,length(new_d)));
h_mpac = zeros(2,length(new_d));
[rx_1_sum,rx_2_sum,rx_real_2_sum] = deal(zeros(length(new_d),length(sig)));

% snr = 20;

%% generate h (channel parameters)
%     x0 = 0.5*lambda;y0 = -0.3*lambda;   %center coordinate of two ants
x0 = 0; y0 = 0;
pos_ant_1 = [-new_d(length(new_d))/2*cos(phi_a), ...
    -new_d(length(new_d))/2*sin(phi_a)]; %(x,y)

S = 2; % number of tx antennas
K = length(phi_sample); % number of probe(position)
CLUSTER = 6; % number of clusters
T = 1e3 + 1; % channel impulse length (unit:ns)
% TAU = 1000; % channel tau(cluster) (unit:ns)
M = 20; % number of subpaths

h_probe = zeros(S,CLUSTER,K,T);

t = linspace(0,1e-6,T);
% phi_sample = linspace(-pi,pi,K);

aoa_main = [65.7 45.6 143.2 32.5 -91.1 -19.2];
aoa_35deg    = [1.5679 4.9447 8.7224 13.0045 17.9492 23.7899 30.9538 40.1824 53.1816 75.4274];      % [1, Table 5.2]
delta_nm_aoa = [aoa_35deg; -aoa_35deg];
delta_nm_aoa = delta_nm_aoa(:);       % these are the same for all users and paths
delta_nm_aoa = repmat(delta_nm_aoa,1,CLUSTER); % a (M x N) matrix
mean_aoas=aoa_main;
aoas=repmat(mean_aoas,M,1)+delta_nm_aoa;
aoas = aoas';

v = 30*1000/(60*60);    % 30km/hour

% this pow and delay are fixed by scenario "micro" maybe should be revised
pow = [-3,-4.3,-5.7,-7.3,-9,-11.4];
pow_mag = sqrt(10.^(pow/10));
delay = [0, 29, 22, 66, 87, 93]; % unit(10ns)

for clu = 1:CLUSTER
    %     for s = 1:S
    PAS_output = generate_PAS_inline(phi_sample,'micro',clu);
    for k = 1:K
        fading_seq= randn(1,T) + 1j*randn(1,T);
        phi_init = -pi + 2*pi*randn(1,M);
        for m = 1:M
            fd(m) = v*cosd(aoas(clu,m) - 120)/lambda;
            % force average power in subpaths
            u_k(m,:) = (1/20) .*exp(1i*(2*pi*fd(m)*t + phi_init(m) ));
        end
        P(clu,k) = pow_mag(clu)*PAS_output.PAS(k);
        h_probe(1,clu,k,:) = P(clu,k)*sum(u_k(m,:),1).*fading_seq;
    end
end
channel_sim = zeros(CLUSTER,K,T);
% h = squeeze( sum( h_probe, 2 ));

%% actually implement channel(h) on spatial correlation
for i = 1:length(new_d)
    %simulate spatial correlation using two antennas with circle
    %         pos_ant_1 = [new_d(i)*cos(phi_a),new_d(i)*sin(phi_a)]; %(x,y)
    pos_ant_2 = pos_ant_1 + [new_d(i)*cos(phi_a),new_d(i)*sin(phi_a)];
    ang_P = zeros(2,length(ideal_phi));
    for j = 1:length(ideal_phi)
        phi = ideal_phi(j) + pi/2;
        %             top = [r*cos(phi+pi/2),r*sin(phi+pi/2)];
        top = [r*cos(phi),r*sin(phi)];
        %             d_1(i,j) = sqrt((pos_ant_1(1) - top(1))^2 + (pos_ant_1(2) - top(2))^2);
        %             d_2(i,j) = sqrt((pos_ant_2(1) - top(1))^2 + (pos_ant_2(2) - top(2))^2);
        [ang_P(1,j),ang_P(2,j)] = generate_ang_of_pattern_v2(top,pos_ant_1,pos_ant_2);
        delta_d(i,j) = d_2(i,j) - d_1(i,j);
        %             h_sig_1(i,j) = exp(1j*(2*pi*(dot(pos_ant_1,top))/lambda)).*sqrt(ideal_PAS(j));
        h_sig_2(i,j) = exp(1j*(2*pi*(dot(pos_ant_2,top))/lambda)).*sqrt(ideal_PAS(j));
        %             rx_1(i,j,:) = P_az_amp(ang_P(1,j))*h_sig_1(i,j)*sig;
        rx_2(i,j,:) = P_az_amp(ang_P(2,j))*h_sig_2(i,j)*sig;
    end
    
    %         rx_1_sum(i,:) = reshape(sum(rx_1(i,:,:),2),1,length(sig));
    rx_2_sum(i,:) = reshape(sum(rx_2(i,:,:),2),1,length(sig));
    
    
    % Corr(:,:,i) = abs(corrcoef(h_new_1(i,:),h_new_2(i,:)));
    % spatial_circle_sig(i) = sum(r1(i,:).*conj(r2(i,:)))/sum(r1(1,:).*conj(r2(1,:)));
    
    %     spatial_circle(i) = sum(h_sig_1(i,:).*conj(h_sig_2(i,:)));
    
    %simulate spatial correlation using two antenna with circle(MPAC),
    %real phi and PAS
    for clu = 1:CLUSTER
        for j = 1:length(phi_sample)
        phi_real = phi_sample(j)+pi/2;
        %%postioner
        % error_top = 0.3*randn*lambda;
        error_top = 0;
        top = [error_top + r*cos(phi_real),error_top + r*sin(phi_real)];   %this errro is used to first positioner
        pos_ant_2 = pos_ant_1 + [new_d(i)*cos(phi_a),new_d(i)*sin(phi_a)];
        
        %error of time alignment on probe
        error = error_para(1) + error_para(2)*randn(1,1);
        % ideal equal pathloss
        % beta angle shows the separate the measurement like SPAC.
        % In MPAC, beta is the different distance of multiple probes to
        % DUT antennas, it rely on the phi_sample (j).
%         h_sig_real_2(i,j) = alpha*exp(1j*(beta(j) + 2*pi*fc*(dot(pos_ant_2,top)/c) )).*sqrt(real_PAS(j));
        [ang_P_1(j),ang_P_2(j)] = generate_ang_of_pattern_v2(top,pos_ant_1,pos_ant_2);
            % h_sig_real_2(i,j) = alpha*exp(1j*(beta(j) + 2*pi*fc*(dot(pos_ant_2,top)/c) )).*sqrt(real_PAS(j));
            % h_probe = zeros(S,CLUSTER,K,T);
            channel_sim(clu,j,:) = h_probe(1,clu,j,:) .* exp(1j*(2*pi*fc*(dot(pos_ant_2,top)/c)));
            % rx_real_2(i,j,:) = P_az_amp(ang_P_2(j))*(h_sig_real_2(i,j)*sig);
            rx_real_cluster_2(i,j,:) = P_az_amp(ang_P_2(j)) .* conv(sig,squeeze(channel_sim(clu,j,:)));
        end
        temp = squeeze(sum(rx_real_cluster_2(i,:,:),2));
        rx_real_cluster(i,clu,:) = [zeros(1,delay(clu)),temp(1:(length(sig) - delay(clu)))'];
    end
    
    %         rx_real_1_sum(i,:) = reshape(sum(rx_real_1(i,:,:),2),1,length(sig));
    rx_real_2_sum(i,:) = squeeze(sum(rx_real_cluster(i,:,:),2));
    
    %randomly choose one point in signal sequence. "100" is chosen randomly
    
    %%num = randi([1,length(sig)],1);
    %         spatial_circle_real_sig_MPAC(i) = sum(rx_real_1(i,:,num).*conj(rx_real_2(i,:,num)));
    
    
    beta1 = 0;
    %using numberical method calculate the spatial correlation
    for j = 1:length(ideal_phi)
        h(i,1,j) = alpha*exp(1j*(beta1)).*sqrt(ideal_PAS(j));
        h(i,2,j) = alpha*exp(1j*(beta1 + 2*pi*(new_d(i)/lambda) * cos((pi/2-ideal_phi(j)) + phi_a) ) ).*sqrt(ideal_PAS(j));
    end
    h1_cal = reshape(h(i,1,:),1,length(ideal_PAS));
    h2_cal = reshape(conj(h(i,2,:)),1,length(ideal_PAS));
    spatial_num(i,:) = sum(h1_cal.*h2_cal);
    
    %     Corr(i,:,:) = abs(corrcoef(squeeze(h(i,1,:)),squeeze(h(i,2,:))));
    
    %using traditional method(ideally equation) calculate
    tau = zeros(1,length(ideal_phi));
    for k = 1:length(ideal_phi)
        tau(k) = d(i)*sin(ideal_phi(k)-phi_a)/c;
        spatial(i) = spatial(i) + exp(-1i*2*pi*fc * tau(k) ).*(ideal_PAS(k));
    end
end

for i = 1:50
    ss(i) = mean(rx_real_2_sum(51-i+1,:).*conj(rx_real_2_sum(51+i-1,:)));
    num = randi([1,length(sig)],1);
    %randomly choose one point in signal sequence. "100" is chosen randomly
    spatial_circle_sig(i) = sum(rx_2(50-i+1,:,num).*conj(rx_2(50+i,:,num)));
    spatial_circle_real_sig_MPAC(i) = sum(rx_real_2(51-i+1,:,num).*conj(rx_real_2(51+i-1,:,num)));
end
res_ss = abs(ss)./abs(ss(1));
spatial = spatial./spatial(1);
spatial_circle_sig = spatial_circle_sig./spatial_circle_sig(1);
spatial_circle_real_sig_MPAC = spatial_circle_real_sig_MPAC./spatial_circle_real_sig_MPAC(1);
% Correaltion = squeeze(Corr(:,1,2));

%calculate the statistical error
stat_MPAC = sqrt(sum(power(abs(spatial_circle_real_sig_MPAC)-abs(spatial_circle_sig),2))/length(spatial_circle_sig));

%     output structure
spatial_output.stat = stat_MPAC;
spatial_output.spatial_sum = res_ss;
spatial_output.spatial_circle_real = spatial_circle_real_sig_MPAC;
spatial_output.spatial_circle = spatial_circle_sig;
spatial_output.spatial_num = spatial_num;
spatial_output.theory = spatial;

% different strength of signal to control probe power

end

% generate PAS data by function:PAS_fun
function PAS_output = generate_PAS_inline(phi_sample,scenario,cluster_index)
switch scenario
    case 'test'
        AOA = -22.5*pi/180;AS = 35*pi/180;
        pow = 0;
        
    case 'micro'
        % % SCME urban micro-cell
        AOA = [0.7,-13.2,146.1,-30.5,-11.4,-1.1] * pi/180;
        AS = 35*pi/180;
        pow = [-3,-4.3,-5.7,-7.3,-9,-11.4];
        
    case 'macro'
        % %SCME urban macro-cell
        AOA = [65.7,45.6,143.2,32.5,-91.1,-19.2]*pi/180;
        AS = 35*pi/180;
        pow = [-3,-5.2,-4.7,-8.2,-12.1,-15.5];
        
    otherwise
        warning('wrong scenario.');
        quit;
end

if strcmp(scenario, 'test')
    cluster_index = 1;
end

amp = 10.^(pow/10);
A = amp./max(amp);
%     PAS_fun = @(x) (1/(2*AS))*sum( A.*exp(-abs(x-AOA)/AS) );
PAS_fun = @(x) (1/(2*AS))* exp(-abs(x - AOA(cluster_index) )/AS);

PAS_output.PAS = zeros(1,length(phi_sample));

for i = 1:length(phi_sample)
    PAS_output.PAS(i) = PAS_fun(phi_sample(i));
end

sumPAS = sum(PAS_output.PAS);
PAS_output.PAS = PAS_output.PAS ./ sumPAS;
PAS_output.pow = pow(cluster_index);
end