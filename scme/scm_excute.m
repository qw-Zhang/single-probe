scmp = scmparset;
scmp.Scenario = 'urban_micro';
% scmp.IntraClusterDsUsed='yes';
scmp.FixedPdpUsed='yes';
scmp.FixedAnglesUsed='yes';
scmp.UniformTimeSampling = 'yes';
% scmp.RandomSeed = 12;
scmp.NumTimeSamples = 1000;
linkp = linkparset;
antp = antparset;

[H, delay, output] = scm(scmp,linkp,antp);
aoa = output.aoas;
aod = output.aods;
phases = output.subpath_phases;
powers = output.path_powers;
%%
max_wavelength = 1;
DIS = linspace(0.001,max_wavelength,100);
corr_para = zeros(1,length(DIS));
for dis_index = 1:length(DIS)
    antp.MsElementPosition = DIS(dis_index);
    [H, delay, output] = scm(scmp,linkp,antp);
    
    H_sum = sum(H,3);
    R1 = squeeze( H(1,1,:) );
    R2 = squeeze( H(2,2,:) );
    corr_temp = corrcoef(R1,R2);
    corr_para(dis_index) = corr_temp(1,2);
end
    
%%
h = H(:,:,:,50);
h1 = reshape(h(1,1,:) + h(2,1,:),1,6);
h2 = reshape(h(1,2,:) + h(2,2,:),1,6);

%%
size(H);

H11 = reshape(H(1,1,:,:),scmp.NumPaths,scmp.NumTimeSamples);

H22 = reshape(H(2,2,:,:),scmp.NumPaths,scmp.NumTimeSamples);

%%
Tx_num = 2; Rx_num = 2;
N = 6; M = 20;
T = linspace(0,1e-4,1e4);
pow_ori = db2mag([-3 -5.2 -4.7 -8.2 -12.1 -15.5]);
pow = pow_ori./sum(pow_ori);
aoa_main = [65.7 45.6 143.2 32.5 -91.1 -19.2];
aod_main = [82 80.5 79.6 98.6 102.1 107.1];
aod_5deg=[0.2236 0.7064 1.2461 1.8578 2.5642 3.3986 4.4220 5.7403 7.5974 10.7753]; % [1, Table 5.2]
        delta_nm_aod = [aod_5deg; -aod_5deg];
        delta_nm_aod=delta_nm_aod(:);   % this (M x 1) vector is the same for all users and paths
        delta_nm_aod=repmat(delta_nm_aod,1,N);  % a (M x N) matrix
        
mean_aods=aod_main;
aods=repmat(mean_aods,M,1)+delta_nm_aod;

aoa_35deg    = [1.5679 4.9447 8.7224 13.0045 17.9492 23.7899 30.9538 40.1824 53.1816 75.4274];      % [1, Table 5.2]
delta_nm_aoa = [aoa_35deg; -aoa_35deg];
delta_nm_aoa = delta_nm_aoa(:);       % these are the same for all users and paths
delta_nm_aoa = repmat(delta_nm_aoa,1,N); % a (M x N) matrix
mean_aoas=aoa_main;
aoas=repmat(mean_aoas,M,1)+delta_nm_aoa;        
        
% delta_aod&aoa are ignored temporally

c = 3e8;
fc = 2.45e9;
lambda = c/fc;

vec_BS = reshape([[0 0], [0 lambda/2]],2,2);
vec_MS = reshape([[0 0], [0 lambda/2]],2,2);

% choose one scenario from "test" "micro" "macro"
scenario = 'micro';
ideal_PAS = generate_PAS(ideal_phi,scenario);
% AS_AOA = 35*180/pi; AOA = aoa_main(n)*180/pi;
% PAS_AOA =@(x) (1/(2*AS_AOA))*(exp(-abs(x-AOA)/AS_AOA));
% 
% AS_AOD = 2*180/pi; AOD = aod_main(n)*180/pi;
% PAS_AOD =@(x) (1/(2*AS_AOD))*(exp(-abs(x-AOD)/AS_AOD));

pattern_BS = @(x) (1);
pattern_MS = @(x) (1);

for u = 1:1
    vec_u = vec_BS(u,:);
    for s = 1:1
        vec_s = vec_MS(s,:);
        for n = 1:1
            for t = 1:length(T)
                for m = 1:M % M = 20
                    AoD = aods(m,n);AoA = aoas(m,n);
                    [x,y] = unit(AoD);
                    vec_AoD = [x,y];
                    [x,y] = unit(AoA);
                    vec_AoA = [x,y];
                    temp = sum( sqrt(pattern_BS(AoD))*sqrt(pattern_MS(AoA)) * ...
                        exp(1i*2*pi/lambda*dot(vec_AoD,(vec_u))) *...
                        exp(1i*2*pi/lambda*dot(vec_AoA,(vec_s))) *...
                        exp(1i*2*pi*dopp(AoA,fc)*t + randn(1) ));
                end
                h(u,s,n,t) = sqrt(pow(n)) * temp;
            end
        end
    end
end



% h_probe = zeros(S,K,CLUSTER,T);
%%
% function PAS_output = generate_PAS(ideal_phi,scenario,cluster_index) 
% generating PAS for every cluster.
fc = 2.45e9;
speed_of_light=2.99792458e8;
wavelength=speed_of_light/fc;
k_CONST = 2*pi/wavelength;
S = 2; % number of tx antennas
K = 16; % number of probe(position)
CLUSTER = 6; % number of clusters
T = 1e4 + 1; % channel impulse length (unit:ns)
TAU = 200; % channel tau(cluster) (unit:5ns)
M = 20; % number of subpaths
MsDirection = 30;
v = 30*1000/(60*60);    % 30km/hour
h_probe = zeros(S,K,CLUSTER,T);
% delta_T = 0.0037;
delta_T = 1e-5;
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
phases= 360*rand(CLUSTER,M);

% according to model and resolution of tau(5ns)
delay = [1,59,42,133,162,186];
%%
max_wavelength = 4;
% DIS = linspace(0,max_wavelength*wavelength,100);
DIS = [0.5];
corr_para = zeros(1,length(DIS));
temp_h = zeros(M,T);
H_test = zeros(2,2,CLUSTER,T);
H_taut = zeros(2,2,TAU,T);
for dis_index = 1:1
    
    phases= 360*rand(CLUSTER,M);
    for u = 1:2 % users number
        %         du = DIS(dis_index) * (u-1) * wavelength;
        du = 0.5 * (u-1) * wavelength;
        for s = 1:2 % tx number
            ds = 0.5 * (s-1) * wavelength;
            for clu = 1:CLUSTER
                for m = 1:M
                    temp_h(m,:) = 1/20 * exp(1j * (k_CONST * ds * sin((aods(clu,m))*pi/180) + (phases(clu,m)*pi/180))) * ...
                        exp(1j * (k_CONST * du * sin((aoas(clu,m))*pi/180))) * ...
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
    
    corr_temp = corrcoef(R1,R2);
    corr_para(dis_index) = corr_temp(1,2);
    % corr_para(dis_index) = mean(R1.*conj(R2));
end

% delay (tau) implement

% x = linspace(0,max_wavelength,length(DIS));
% norm_corr = corr_para ./ corr_para(1);
% figure;plot(x,abs(norm_corr));
% axis([0 max_wavelength 0 1]);

%%
% function PAS_output = generate_PAS(ideal_phi,scenario,cluster_index) 
% generating PAS for every cluster.

fc = 2.45e9;
S = 2; % number of tx antennas
K = 16; % number of probe(position)
CLUSTER = 6; % number of clusters
T = 1e3 + 1; % channel impulse length (unit:ns)
TAU = 1000; % channel tau(cluster) (unit:ns)
M = 20; % number of subpaths

h_probe = zeros(S,K,CLUSTER,T);

t = linspace(0,1e-6,T);
phi_sample = linspace(-pi,pi,K);

aoa_main = [65.7 45.6 143.2 32.5 -91.1 -19.2];
aoa_35deg    = [1.5679 4.9447 8.7224 13.0045 17.9492 23.7899 30.9538 40.1824 53.1816 75.4274];      % [1, Table 5.2]
delta_nm_aoa = [aoa_35deg; -aoa_35deg];
delta_nm_aoa = delta_nm_aoa(:);       % these are the same for all users and paths
delta_nm_aoa = repmat(delta_nm_aoa,1,CLUSTER); % a (M x N) matrix
mean_aoas=aoa_main;
aoas=repmat(mean_aoas,M,1)+delta_nm_aoa;
aoas = aoas';

pow = [-3,-4.3,-5.7,-7.3,-9,-11.4];
pow_mag = sqrt(10.^(pow/10));

for clu = 1:CLUSTER
%     for s = 1:S
        PAS_output = generate_PAS(phi_sample,'micro',clu);
        for k = 1:K
            fading_seq= rayleigh(T);
            phi_init = -pi + 2*pi*randn(1,M);
            for m = 1:M
                fd(m) = dopp(120,aoas(clu,m),fc);
                % force average power in subpaths
                u_k(m,:) = (1/20) .*exp(1i*(2*pi*fd(m)*t + phi_init(m) ));
            end
            P(clu,k) = pow_mag(clu)*PAS_output.PAS(k);
            h_probe(1,k,clu,:) = P(clu,k)*sum(u_k(m,:),1).*fading_seq;
        end
end

h = squeeze( sum( h_probe, 2 ));

% implement spatial character


% functions
function [x,y] = unit(ang)
    x = cosd(ang);
    y = sind(ang);
end

function fading_seq = rayleigh(len_channel)
    fading_seq = randn(1,len_channel) + 1j*randn(1,len_channel);
end

function fd = dopp(ang_mov,ang_aoa,fc)
    c = 3e8;
    lambda = c/fc;
    v = 30*1000/(60*60);    % 30km/hour
    fd = v*cosd(ang_aoa - ang_mov)/lambda;
end

% generate PAS data by function:PAS_fun
function PAS_output = generate_PAS(phi_sample,scenario,cluster_index)
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