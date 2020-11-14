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
size(H)

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
t = 0:10e-9:1e-6;
tau = delay;
phase_random = -pi + 2*pi*rand(6,20,36);
c = 3e8;
fre = 2e9;
lambda = c / fre;
doppler_fre = cos( deg2rad( reshape(aoa,6,20) )) * 10 / lambda; %   speed = 10
temp = 0;
power_weight = ones(36,6);
for u = 1:2
    ant_gain = ones(1,360);
    for k = 1:36
        for n = 1:6
            for m = 1:20
                temp = temp + exp(1j*2*pi*doppler_fre(n,m) + phase_random(n,m,k));
            end
            hota_ukn(u,k,n) = sqrt(powers(n)) * sqrt(ant_gain(round(aod(1,n,m)))) *temp * sqrt(power_weight(k,n));
        end
    end
end
% delayseq slove 'n'