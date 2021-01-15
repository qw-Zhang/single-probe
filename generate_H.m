% h_probe = zeros(S,K,CLUSTER,T);
%%
% function PAS_output = generate_PAS(ideal_phi,scenario,cluster_index) 
% generating PAS for every cluster.
function [H,PAS,output_para] = generate_H(phi_sample,real_PAS)
fc = 2.45e9;
speed_of_light=2.99792458e8;
wavelength=speed_of_light/fc;
k_CONST = 2*pi/wavelength;
S = 2; % number of tx antennas
K = 16; % number of probe(position)
CLUSTER = 6; % number of clusters
T = 1e3 + 1; % channel impulse length (unit:ns)
TAU = 200; % channel tau(cluster) (unit:5ns)
M = 20; % number of subpaths
MsDirection = 30;
v = 30*1000/(60*60);    % 30km/hour
h_probe = zeros(S,K,CLUSTER,T);
delta_T = 0.0037;
% delta_T = 2e-3;
t = linspace(0,T*delta_T,T);
phi_sample = linspace(-pi,pi,K);

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

% if strcmpi(IntraClusterDsUsed,'yes')
%     NumSubPathsPerMidpath=[6,6,4,4]; MidPathOrder=
%     [1,2,3,4,19,20,5,6,7,8,17,18,9,10,15,16,11,12,13,14];
%     MidPathDelays=[0  5.8e-009  1.35e-008  2.76e-008];    % relative to
%     path/cluster delay MidPathPowers=[6/20 6/20 4/20 4/20];  % relative
%     to path power
% end

pow = [-3,-5.7,-4.3,-7.3,-9,-11.4];
pow_mag = sqrt(10.^(pow/10));
pow_mag_norm = pow_mag ./ sum(pow_mag);

% according to model and resolution of tau(5ns)
delay = [1,59,42,133,162,186];

max_wavelength = 4;
DIS = linspace(0,max_wavelength*wavelength,20);
% DIS = [0.5];
corr_res = zeros(1,length(DIS));
temp_h = zeros(M,T);
H_test = zeros(2,2,CLUSTER,T);
H_taut = zeros(2,2,TAU,T);
phases= 360*rand(CLUSTER,M);

% just for test
ideal_phi = linspace(-pi,pi,3600);
% real_phi as a parameter
phi_sample = linspace(-pi,pi,16);
real_PAS = zeros(1,length(phi_sample));
rate = floor(length(ideal_phi)/length(phi_sample));
% rate = length(ideal_phi)/length(real_phi);

%scenario parameters
% phi_a = 0*pi/180;
scenario = 'micro';
ideal_PAS = generate_PAS1(ideal_phi,scenario);

% real probe scenario
for i = 1:length(phi_sample)
    real_PAS(i) = sum(ideal_PAS(rate*(i-1)+1:rate*i))/rate;
    phi_sample(i) = sum(ideal_phi(rate*(i-1)+1:rate*i))/rate;
end
phi_sample_deg = phi_sample * 180 / pi;
real_sum_PAS = sum(real_PAS);
real_PAS = real_PAS/real_sum_PAS;


% for dis_index = 1:length(DIS)
for dis_index = 10:10
    for u = 1:2 % users number
          du = DIS(dis_index) * (u-1) * wavelength;
%         du = 0.5 * (u-1) * wavelength;
        for s = 1:2 % tx number
            ds = 0.5 * (s-1) * wavelength;
            for clu = 1:CLUSTER
                for m = 1:length(phi_sample)
                    temp_h(m,:) = real_PAS(m) * exp(1j * (k_CONST * ds * sin((aods(clu,m))*pi/180) + (phases(clu,m)*pi/180))) * ...
                                         exp(1j * (k_CONST * sin(phi_sample_deg(m)*pi/180))) * ...
                                         exp(1j * k_CONST * v * cos((aoas(clu,m) - MsDirection)*pi/180) .* t);
                end
                H_test(u,s,clu,:) = pow_mag_norm(clu) .* sum(temp_h,1);
                
                for t_index = 1:T
                    H_taut(u,s,delay(clu),t_index) = H_test(u,s,clu,t_index);
                end
            end
        end
    end
    
    R1_temp = squeeze( H_taut(1,1,:,:) + H_taut(2,1,:,:) );
%     R1_t = 
    R2_temp = squeeze( H_taut(1,2,:,:) + H_taut(2,2,:,:) );
    
%     for t_index = 1:T
%        rxSig1_interp(t_index,:) = conv(txSig_inter,R1_temp(:,t_index));
%        rxSig1_temp(t_index,:) = resample(rxSig1_interp(t_index,1:size(txSig_inter,1)),1,10);
%        rxSig2_interp(t_index,:) = conv(txSig_inter,R2_temp(:,t_index));
%        rxSig2_temp(t_index,:) = resample(rxSig2_interp(t_index,1:size(txSig_inter,1)),1,10);
%     end
%     
%     rxSig1_res = reshape(rxSig1_temp, size(rxSig1_temp,1) * size(rxSig1_temp,2),1);
%     rxSig2_res = reshape(rxSig2_temp, size(rxSig2_temp,1) * size(rxSig2_temp,2),1);
% 
%     corr_temp = corrcoef(rxSig1_res(1:74000,1),rxSig2_res(1:74000,1));
%     corr_res(dis_index) = corr_temp(1,2);
    
end

% delay (tau) implement

% x = linspace(0,max_wavelength,length(DIS));
% norm_corr = corr_para ./ corr_para(1);
% figure;plot(x,abs(norm_corr));
% axis([0 max_wavelength 0 1]);