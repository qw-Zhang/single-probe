% h_probe = zeros(S,K,CLUSTER,T);
%%
% function PAS_output = generate_PAS(ideal_phi,scenario,cluster_index)
% generating PAS for every cluster.
function [h,output_para] = generate_H(phi_sample,scenario)
fc = 2.45e9;
speed_of_light=2.99792458e8;
wavelength=speed_of_light/fc;
k_CONST = 2*pi/wavelength;
S = 2; % number of tx antennas
K = length(phi_sample); % number of probe(position)
CLUSTER = 6; % number of clusters
T = 5e2 + 1; % channel impulse length (unit:ns)
% TAU = 200; % channel tau(cluster) (unit:5ns)
M = 20; % number of subpaths
MsDirection = 60;
v = 30*1000/(60*60);    % 30km/hour
h_probe = zeros(S,K,CLUSTER,T);
delta_T = 0.0037;
% delta_T = 2e-3;
t = linspace(0,T*delta_T,T);

switch scenario
    case 'micro'
        aod_5deg=[0.2236 0.7064 1.2461 1.8578 2.5642 3.3986 4.4220 5.7403 7.5974 10.7753]; % [1, Table 5.2]
        delta_nm_aod = [aod_5deg; -aod_5deg];
        delta_nm_aod=delta_nm_aod(:);   % this (M x 1) vector is the same for all users and paths
        delta_nm_aod=repmat(delta_nm_aod,1,CLUSTER);  % a (M x N) matrix
        mean_aods=[6.6100 50.8297 14.1360 38.3972 6.6690 40.2849];
        aods=repmat(mean_aods,M,1)+delta_nm_aod;
        aods = aods';
        
        aoa_main = [0.6966 146.0669 -13.2268  -30.5485 -11.4412 -1.0587];
        aoa_35deg    = [1.5679 4.9447 8.7224 13.0045 17.9492 23.7899 30.9538 40.1824 53.1816 75.4274];      % [1, Table 5.2]
        delta_nm_aoa = [aoa_35deg; -aoa_35deg];
        delta_nm_aoa = delta_nm_aoa(:);       % these are the same for all users and paths
        delta_nm_aoa = repmat(delta_nm_aoa,1,CLUSTER); % a (M x N) matrix
        mean_aoas=aoa_main;
        aoas=repmat(mean_aoas,M,1)+delta_nm_aoa;
        aoas = aoas';
        
        pow = [-3,-5.7,-4.3,-7.3,-9,-11.4];
        delay = [1,30,22,67,82,94];
        %         delay = [1,80,180,300,400,800];
        
    case 'macro'
        aod_2deg     = [0.0894 0.2826 0.4984 0.7431 1.0257 1.3594 1.7688 2.2961 3.0389 4.3101];     % [1, Table 5.2]
        delta_nm_aod = [aod_2deg; -aod_2deg];
        delta_nm_aod = delta_nm_aod(:);   % this (M x 1) vector is the same for all users and paths
        delta_nm_aod = repmat(delta_nm_aod,1,CLUSTER);  % a (M x N) matrix
        mean_aods=[-127.2788  -136.8071 -129.9678 -96.2155 -159.5999 173.1860];
        aods=repmat(mean_aods,M,1)+delta_nm_aod;
        aods = aods';
        
        aoa_35deg    = [1.5679 4.9447 8.7224 13.0045 17.9492 23.7899 30.9538 40.1824 53.1816 75.4274];      % [1, Table 5.2]
        delta_nm_aoa = [aoa_35deg; -aoa_35deg];
        delta_nm_aoa = delta_nm_aoa(:);       % these are the same for all users and paths
        delta_nm_aoa = repmat(delta_nm_aoa,1,CLUSTER); % a (M x N) matrix
        mean_aoas=[65.7489 143.1863 45.6454 32.5131 -91.0551 -19.1657];
        aoas=repmat(mean_aoas,M,1)+delta_nm_aoa;
        aoas = aoas';
        
        taus_sorted=[0  2.527e-007  3.600e-007  1.0387e-006  2.7300e-006  4.5977e-006];
        pow = [-3,-5.2,-4.7,-8.2,-12.1,-15.5];
        delay = [1,37,27,105,274,461];
        
    case 'test'
        CLUSTER = 1;
        
        aod_2deg     = [0.0894 0.2826 0.4984 0.7431 1.0257 1.3594 1.7688 2.2961 3.0389 4.3101];     % [1, Table 5.2]
        delta_nm_aod = [aod_2deg; -aod_2deg];
        delta_nm_aod = delta_nm_aod(:);   % this (M x 1) vector is the same for all users and paths
        delta_nm_aod = repmat(delta_nm_aod,1,CLUSTER);  % a (M x N) matrix
        mean_aods=[-127.2788];
        aods=repmat(mean_aods,M,1)+delta_nm_aod;
        aods = aods';
        
        
        AOA = -22.5*pi/180;AS = 35*pi/180;
        aoa_35deg    = [1.5679 4.9447 8.7224 13.0045 17.9492 23.7899 30.9538 40.1824 53.1816 75.4274];      % [1, Table 5.2]
        delta_nm_aoa = [aoa_35deg; -aoa_35deg];
        delta_nm_aoa = delta_nm_aoa(:);       % these are the same for all users and paths
        delta_nm_aoa = repmat(delta_nm_aoa,1,CLUSTER); % a (M x N) matrix
        mean_aoas=[-22.5];
        aoas=repmat(mean_aoas,M,1)+delta_nm_aoa;
        aoas = aoas';
        pow = 0;
        delay = [1];
end

P = zeros(CLUSTER,K);
% according to model and resolution of tau(10ns)

subpath = zeros(M,T);
h_test = zeros(2,CLUSTER,K,T);


phi_sample_deg = phi_sample .* 180 / pi;


for s = 1:2 % tx number
    ds = 0.5 * (s-1) * wavelength;
    for clu = 1:CLUSTER
        %         h_temp = (randn(K,length(t)) + 1i*randn(K,length(t)));
        PAS_output_clu = generate_PAS_inline(phi_sample,scenario,clu);
        %             phases = 360*rand(K,M);
        phases = -180 + 360*rand(K,M);
        for k = 1:K
            for m = 1:M
                subpath(m,:) = 1/20 * exp(1j * (k_CONST * v * cos((aoas(clu,m) - MsDirection)*pi/180) .* t + (phases(k,m)*pi/180)) );
            end
            amp_mag = db2mag(PAS_output_clu.pow);
            pow_mag_norm = power( (amp_mag ./ db2mag(PAS_output_clu.pow_init)),2);
            P(clu,k) = PAS_output_clu.PAS(k) * pow_mag_norm;
            %             P_root(clu,k) = PAS_output_clu.PAS(k) * nthroot(pow_mag_norm,4);
            h_test(s,clu,k,:) = sqrt(P(clu,k)) .* sum(subpath,1);
            % h_test(s,clu,k,:) = P(clu,k) .* h_temp(k,:);
        end
    end
end

% return sum of all clusters
PAS_output_temp = repmat(PAS_output_clu, CLUSTER, 1 );
ideal_length = length(PAS_output_clu.ideal_PAS);
P_ideal_sum_temp = zeros(CLUSTER,ideal_length);
P_real_sum_temp = zeros(CLUSTER, length(PAS_output_clu.PAS));
for clu = 1:CLUSTER
    PAS_output_temp(clu) = generate_PAS_inline(phi_sample,scenario,clu);
    %                 amp_mag = sqrt(10.^(PAS_output_clu.pow/10));
    amp_mag = db2mag(PAS_output_temp(clu).pow);
    pow_mag_norm = power( (amp_mag ./ db2mag(PAS_output_temp(clu).pow_init)),2);
    P_ideal_sum_temp(clu,:) = PAS_output_temp(clu).ideal_PAS * pow_mag_norm;
%     P_real_sum_temp(clu,:) = PAS_output_clu.PAS * pow_mag_norm;
	P_real_sum_temp(clu,:) = PAS_output_temp(clu).PAS .* pow_mag_norm;
end
P_ideal_sum = sum(P_ideal_sum_temp,1);
P_real_sum = sum(P_real_sum_temp,1);

h.h_ori = h_test;
% h.h1 = squeeze(sum(h_test,1));
% h.h2 = squeeze(sum(h_test,1));

% h.h1 = reshape(sum(h_test,1),CLUSTER,K,T);
h.h1 = reshape(h_test(1,:,:,:),CLUSTER,K,T);

% h.h1 = reshape(sum(h_test,1),CLUSTER,K,T);
output_para.real_PAS_sum = P_real_sum;
output_para.ideal_PAS_sum = P_ideal_sum;
output_para.pow = P;
output_para.h_aoas = aoas;
output_para.h_aods = aods;
output_para.delay = delay;
output_para.delta_T = delta_T;
end

%%
% generate PAS data by function:PAS_fun
function PAS_output = generate_PAS_inline(phi_sample,scenario,cluster_index)
switch scenario
    case 'test'
        AOA = -22.5*pi/180;AS = 35*pi/180;
        pow = 0;
        
    case 'micro'
        % % SCME urban micro-cell
        AOA = [0.7,-13.2,146.1,-30.5,-11.4,-1.1] * pi/180;
%         AOA = [0.7,-13.2,146.1,-30.5,-21.4,30.1] * pi/180;
        AS = 35*pi/180;
        pow_table = [ [-3 -5.2 -7]; [-4.3 -6.5 -8.3]; [-5.7 -7.9 -9.7]; [-7.3 -9.5 -11.3]; [-9 -11.2 -13]; [-11.4 -13.6 -15.4] ];
        pow = mean(pow_table,2);
        
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

% PAS_fun = @(x) (1/(2*AS))* exp(-abs(x - AOA(cluster_index) )/AS);
PAS_fun = @(x) (1/(2*AS))* exp(-abs(x - 0 )/AS);

ideal_phi = linspace(-pi,pi,3600);
ideal_PAS_temp = zeros(size(ideal_phi));
% ideal_PAS = generate_PAS1(ideal_phi,scenario);
%
% % real probe scenario
% for i = 1:length(phi_sample)
%     real_PAS(i) = sum(ideal_PAS(rate*(i-1)+1:rate*i))/rate;
%     phi_sample(i) = sum(ideal_phi(rate*(i-1)+1:rate*i))/rate;
% end
% real_sum_PAS = sum(real_PAS);
% real_PAS = real_PAS/real_sum_PAS;

for ideal_index = 1:length(ideal_phi)
    ideal_PAS_temp(ideal_index) = PAS_fun( ideal_phi(ideal_index) );
end

move_step = round(AOA(cluster_index) * 180 / pi * 10);
ideal_PAS = circshift(ideal_PAS_temp,move_step);
ideal_PAS = ideal_PAS ./ max(ideal_PAS);
rate = floor(length(ideal_phi)/length(phi_sample));
% real probe scenario
% for i = 1:length(phi_sample)
%     PAS_output.PAS(i) = sum(ideal_PAS(rate*(i-1)+1:rate*i))/rate;
%     phi_sample(i) = sum(ideal_phi(rate*(i-1)+1:rate*i))/rate;
% end
real_PAS_temp = downsample(ideal_PAS,rate);
PAS_output.phi_sample = downsample(ideal_phi,rate);
real_sum_PAS = sum(real_PAS_temp);
% real_PAS = real_PAS_temp./real_sum_PAS;
real_PAS = real_PAS_temp;

% amp = 10.^(pow/10);
% A = amp./max(amp);
%     PAS_fun = @(x) (1/(2*AS))*sum( A.*exp(-abs(x-AOA)/AS) );


% PAS_output.PAS = zeros(1,length(phi_sample));
%
% for i = 1:length(phi_sample)
%     PAS_output.PAS(i) = PAS_fun(phi_sample(i));
% end

% just for test
PAS_output.ideal_PAS = ideal_PAS;
PAS_output.pow_init = pow(1);
PAS_output.PAS = real_PAS;
PAS_output.pow = pow(cluster_index);
end