scmp = scmparset;
scmp.RandomSeed = 0;
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
% normal mimo  (u:element of BS,s:element of MS,n:path,m:subpath of path,():function,[]:value)
% h_u,s,n(t,tau) = sqrt(p[n])*sum_m( sqrt(antenna_gain( aod[n,m] )) * sqrt(antenna_gain( aoa[n,m] ))
%                                   * exp(1j*2*pi* doppler_fre[n,m] *t + phase_random[n,m]) * dirac( tau[n] ) )
%
% SPAC scheme     (u:element of BS,k:position of probe,n:path,m:subpath of path,():function,[]:value)
% h_ota_u,k,n(t,tau) = sqrt(p[n])*sum_m( sqrt(antenna_gain( aod[n,m] ))
%                                   * exp(1j*2*pi* doppler_fre[n,m] *t + phase_random[n,m,k]) )
%                                   * dirac( tau[n] ) * sqrt(power_weight[k,n])
% h_k(t,tau) = sum_u( sum_n( h_ota_u,k,n ) )

ideal_phi = linspace(-pi,pi,3600);
AOA = [0.7,-13.2,146.1,-30.5,-11.4,-1.1]*pi/180;AS = 35*pi/180;
step_mov = floor([0.7,-13.2,146.1,-30.5,-11.4,-1.1]);
pow = [-3,-4.3,-5.7,-7.3,-9,-11.4];

%ideal scenario
amp = 10.^(pow/10);
A = amp./max(amp);
PAS_norm = zeros(1,length(ideal_phi));
for i = 1:length(ideal_phi)
    % Lap distr
    PAS_norm(i) = 1*(1/(2*AS))*(exp(-abs(ideal_phi(i)-0)/AS));
end

PAS_move_temp = zeros(length(step_mov),length(PAS_norm));
for j = 1:length(step_mov)
    PAS_norm = A(j)*PAS_norm;
    if(step_mov(j) < 0)
        PAS_move_temp(j,:) = [PAS_norm(abs(step_mov(j))+1:length(PAS_norm)),PAS_norm(1:abs(step_mov(j)))];
    else
        PAS_move_temp(j,:) = [PAS_norm(length(PAS_norm) - abs(step_mov(j)) + 1:length(PAS_norm))...
            ,PAS_norm(1:length(PAS_norm) - abs(step_mov(j)))];
    end
end
PAS = sum(PAS_move_temp,1);
%using for normalized PAS
sum_PAS = sum(PAS);
PAS = PAS/sum_PAS;
phi_sample = ones(1,36);
rate = length(ideal_phi) / length(phi_sample);
for i = 1:length(phi_sample)
    real_PAS(i) = sum(PAS(rate*(i-1)+1:rate*i))/rate;
    phi_sample(i) = sum(ideal_phi(rate*(i-1)+1:rate*i))/rate;
end
real_sum_PAS = sum(real_PAS);
real_PAS = real_PAS/real_sum_PAS;


time = 0:1e-8:1e-6;
tau = round(delay*1e8);
phase_random = -pi + 2*pi*rand(6,20,36);
c = 3e8;
fre = 2e9;
lambda = c / fre;
doppler_fre = cos(deg2rad( reshape(aoa,6,20) )) * 10 / lambda; %   speed = 10
temp = 0;
power_weight = ones(36,6);
hota_ukn = zeros(length(time),6,2,36);
t = 2;
% for t = 1:length(time)
    for u = 1:2
        ant_gain = ones(1,360);
        for pos = 1:36
            for path = 1:6
                for m = 1:20
                    temp = temp + sqrt(ant_gain(180+round(aod(1,path,m)))) * exp(1i*(2*pi*doppler_fre(path,m)*time(t) + phase_random(path,m,pos)));
                end
                hota_ukn(t,path,u,pos) = db2mag(-3.04)*sqrt(powers(path)) * temp * sqrt(power_weight(pos,path))./20;
            end
        end
    end
% end
channel = zeros(36,6,1000); %  channel matrix

hh = reshape(hota_ukn(2,:,:,:),6,2,36);
hh1 = reshape(sum(hh,2),6,36);
hh2 = hh1.*sqrt(real_PAS);


for pos = 1:36
for tt = 1:length(tau)
    delayseq
    channel(pos,tt,:) = 
% delayseq slove 'n'