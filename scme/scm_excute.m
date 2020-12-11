scmp = scmparset;
scmp.Scenario = 'urban_micro';
% scmp.IntraClusterDsUsed='yes';
scmp.FixedPdpUsed='yes';
scmp.FixedAnglesUsed='yes';
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
pow = db2mag([-3 -5.2 -4.7 -8.2 -12.1 15.5]);
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



function [x,y] = unit(ang)
    x = cos(deg2rad(ang));
    y = sin(deg2rad(ang));
end

function fd = dopp(ang_v,fc)
    c = 3e8;
    lambda = c/fc;
    v = 30*1000/(60*60);    % 30km/hour
    fd = v*cos(deg2rad(ang_v-120))/lambda;
end