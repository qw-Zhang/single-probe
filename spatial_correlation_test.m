%%
%generate Laplace distribution
ideal_phi = linspace(0,2*pi,100);
miu = [0.7,-13.2,146.1,-30.5,-11.4,-1.1]*pi/180;b = 35*pi/180;
figure;hold on;
for i = 3:3
    Lap(i,:) =(1/(2*b))*exp(-abs(ideal_phi-miu(i))/b);
    % integral(Lap,-pi,pi)
    % plot(x,Lap(i,:));
end
pow = [-3,-4.3,-5.7,-7.3,-9,-11.4];
amp = 10.^(pow/10);
amp = amp./max(amp);
for i = 3:3
    LAP(i,:) = amp(i).*Lap(i,:);
    plot(ideal_phi,LAP(i,:));
end
Lap_fin = sum(LAP); Lap_fin = Lap_fin/sum(Lap_fin);
%% 
fc = 2.535e8;
c = 3e8;
lambda = c/fc;
d = linspace(0,lambda,100);

phi_a = 0*pi/180;
scenario = 'micro';

if isequal(scenario,'test')
    AOA = -22.5*pi/180;AS = 35*pi/180;
    PAS =@(x) (1/(2*AS))*(exp(-abs(x-AOA)/AS));
end
if isequal(scenario,'micro')
    % % SCME urban micro-cell
    AOA = [0.7,-13.2,146.1,-30.5,-11.4,-1.1]*pi/180;AS = 35*pi/180;phi_a = 0*pi/180;
    pow = [-3,-4.3,-5.7,-7.3,-9,-11.4];
    amp = 10.^(pow/10);
    A = amp./max(amp);
    PAS =@(x) (1/(2*AS))*(A(1)*exp(-abs(x-AOA(1))/AS)+A(2)*exp(-abs(x-AOA(2))/AS)+A(3)*exp(-abs(x-AOA(3))/AS)+...
        A(4)*exp(-abs(x-AOA(4))/AS)+A(5)*exp(-abs(x-AOA(5))/AS)+A(6)*exp(-abs(x-AOA(6))/AS));
elseif isequal(scenario,'macro')
    % %SCME urban macro-cell
    AOA = [65.7,45.6,143.2,32.5,-91.1,-19.2]*pi/180;AS = 35*pi/180;phi_a = 0*pi/180;
    pow = [-3,-5.2,-4.7,-8.2,-12.1,-15.5];
    amp = 10.^(pow/10);
    A = amp./max(amp);
    PAS =@(x) (1/(2*AS))*(A(1)*exp(-abs(x-AOA(1))/AS)+A(2)*exp(-abs(x-AOA(2))/AS)+A(3)*exp(-abs(x-AOA(3))/AS)+...
        A(4)*exp(-abs(x-AOA(4))/AS)+A(5)*exp(-abs(x-AOA(5))/AS)+A(6)*exp(-abs(x-AOA(6))/AS));
end

spatial = zeros(1,length(d));
sum_PAS = integral(PAS,-pi,pi);
for i = 1:100
fun_int = @(phi_p) exp(-1i*2*pi*(d(i)/lambda)*sin(phi_p - phi_a)).*PAS(phi_p)/sum_PAS;
spatial(i) = integral(fun_int,-pi,pi);
end

figure;plot(d/lambda,abs(spatial));
% grid on;

%%
%PAS continous->distrete
fc = 2.535e9;
c = 3e8;
lambda = c/fc;
d = linspace(0,lambda,100);
ideal_phi = linspace(-pi,pi,180);
PAS = zeros(1,length(ideal_phi));

real_phi = linspace(-pi,pi,18);
real_PAS = zeros(1,length(real_phi));

rate = length(ideal_phi)/length(real_phi);

phi_a = 0*pi/180;
scenario = 'micro';

switch scenario
    case 'test'
        AOA = -22.5*pi/180;AS = 35*pi/180;
        pow = 0;
        
    case 'micro'
        % % SCME urban micro-cell
        AOA = [0.7,-13.2,146.1,-30.5,-11.4,-1.1]*pi/180;AS = 35*pi/180;
        pow = [-3,-4.3,-5.7,-7.3,-9,-11.4];
        
    case 'macro'
        % %SCME urban macro-cell
        AOA = [65.7,45.6,143.2,32.5,-91.1,-19.2]*pi/180;AS = 35*pi/180;
        pow = [-3,-5.2,-4.7,-8.2,-12.1,-15.5];
        
    otherwise
        error('wrong scenario.');
end

%ideal scenario
amp = 10.^(pow/10);
A = amp./max(amp);
for i = 1:length(ideal_phi)
    for j = 1:length(AOA)
        PAS(i) = PAS(i) + A(j)*(1/(2*AS))*(exp(-abs(ideal_phi(i)-AOA(j))/AS));
    end
end

%real probe scenario
for i = 1:length(real_phi)
%     for j = 1:length(AOA)
%         real_PAS(i) = real_PAS(i) + A(j)*(1/(2*AS))*(exp(-abs(real_phi(i)-AOA(j))/AS));
%     end
    
    real_PAS(i) = sum(PAS(rate*(i-1)+1:rate*i))/rate;
    real_phi(i) = sum(ideal_phi(rate*(i-1)+1:rate*i))/rate;
end

spatial = zeros(2,length(d));
% sum_PAS = integral(PAS,-pi,pi);
sum_PAS = sum(PAS);
real_sum_PAS = sum(real_PAS);
error_tk = 1e-9*randn(1,length(real_phi));
% error_tk = zeros(1,length(real_phi));
for i = 1:length(d)
    for j = 1:length(ideal_phi)
        spatial(1,i) = spatial(1,i) + exp(-1i*2*pi*(d(i)/lambda)*sin(ideal_phi(j) - phi_a)).*PAS(j)/sum_PAS;
    end
    
    for k = 1:length(real_phi)
        tau(k) = d(i)*sin(real_phi(k)-phi_a)/c;
        spatial(2,i) = spatial(2,i) + exp(-1i*2*pi*fc * (tau(k) + error_tk(k)) ).*real_PAS(k)/real_sum_PAS;
    end
end

% stat = sqrt(sum((abs(spatial(2,:)-spatial(1,:)).^2))/length(spatial(1,:)))
stat = sum(abs(spatial(2,:)-spatial(1,:)))/length(spatial(1,:))

figure;plot(d/lambda,abs(spatial(1,:)),'b');
hold on;
plot(d/lambda,abs(spatial(2,:)),'r');
axis([0 1 0 1]);
xlabel('Antenna Separation in wavelength');
ylabel('Spatial Correlation')
grid on;

%%
%2020.01.07
%test with h1 h2 without the euqalation of saptial correlation
%using three method generate spatial correlation:
%   1. traditional formula method
%   2. numberical calculation method
%   3. simulation two antennas with circle method
clear;
fc = 2.535e9;
c = 3e8;
lambda = c/fc;
d = linspace(0,lambda,100);
new_d = linspace(0,lambda/2,100);
ideal_phi = linspace(-pi,pi,180);
PAS = zeros(1,length(ideal_phi));

real_phi = linspace(-pi,pi,18);
real_PAS = zeros(1,length(real_phi));

rate = length(ideal_phi)/length(real_phi);

phi_a = 60*pi/180;
scenario = 'micro';

switch scenario
    case 'test'
        AOA = -22.5*pi/180;AS = 35*pi/180;
        pow = 0;
        
    case 'micro'
        % % SCME urban micro-cell
        AOA = [0.7,-13.2,146.1,-30.5,-11.4,-1.1]*pi/180;AS = 35*pi/180;
        pow = [-3,-4.3,-5.7,-7.3,-9,-11.4];
        
    case 'macro'
        % %SCME urban macro-cell
        AOA = [65.7,45.6,143.2,32.5,-91.1,-19.2]*pi/180;AS = 35*pi/180;
        pow = [-3,-5.2,-4.7,-8.2,-12.1,-15.5];
        
    otherwise
        error('wrong scenario.');
end

%ideal scenario
amp = 10.^(pow/10);
A = amp./max(amp);
for i = 1:length(ideal_phi)
    for j = 1:length(AOA)
        PAS(i) = PAS(i) + A(j)*(1/(2*AS))*(exp(-abs(ideal_phi(i)-AOA(j))/AS));
    end
end

%real probe scenario
for i = 1:length(real_phi)
    real_PAS(i) = sum(PAS(rate*(i-1)+1:rate*i))/rate;
    real_phi(i) = sum(ideal_phi(rate*(i-1)+1:rate*i))/rate;
end

spatial = zeros(2,length(d));
% beta = -pi + (pi--pi).*rand(1,length(ideal_phi));
% beta = -pi + (pi--pi).*rand(1);
beta = 0;
alpha = 1;

%using for normalized PAS
sum_PAS = sum(PAS);
real_sum_PAS = sum(real_PAS);

%random error
% error_tk = 1e-11*randn(1,length(real_phi));
error_tk = zeros(1,length(real_phi));


%using numberical method calculate the spatial correlation
spatial_num = zeros(length(d),1);

%simulate radius of circle
r = 1;

for i = 1:length(d)
    %simulate spatial correlation using two antennas with circle 
    for j = 1:length(ideal_phi)
        phi_1 = (ideal_phi(j) - phi_a);
        phi_2 = pi - phi_1;
        d_1(j) = sqrt( r^2 + new_d(i)^2 - 2*new_d(i)*r*cos( phi_1 - pi/2 ) );
        d_2(j) = sqrt( r^2 + new_d(i)^2 - 2*new_d(i)*r*cos( phi_2 + pi/2 ) );
        delta_l(j) = d_2(j) - d_1(j);
        h_new_1(i,j) = alpha*exp(1j*(beta)).*sqrt(PAS(j)/sum_PAS);
        h_new_2(i,j) = alpha*exp(1j*(beta + 2*pi*delta_l(j)/lambda)).*sqrt(PAS(j)/sum_PAS);
    end
    spatial_circle(i) = sum(h_new_1(i,:).*conj(h_new_2(i,:)));
       
    %using numberical method calculate the spatial correlation
    for j = 1:length(ideal_phi)
        h(i,1,j) = alpha*exp(1j*(beta)).*sqrt(PAS(j)/sum_PAS);
        h(i,2,j) = alpha*exp(1j*(beta + 2*pi*(2*new_d(i)/lambda) * cos((pi/2-ideal_phi(j)) + phi_a) ) ).*sqrt(PAS(j)/sum_PAS);
    end
    h1_cal = reshape(h(i,1,:),1,length(PAS));
    h2_cal = reshape(conj(h(i,2,:)),1,length(PAS));
    spatial_num(i,:) = sum(h1_cal.*h2_cal);
    
    %using traditional method(ideally equation) calculate
    for k = 1:length(real_phi)
        tau(k) = d(i)*sin(real_phi(k)-phi_a)/c;
        spatial(2,i) = spatial(2,i) + exp(-1i*2*pi*fc * (tau(k) + error_tk(k)) ).*real_PAS(k)/real_sum_PAS;
    end
end

%calculate the statistical error
% stat = sum(abs(spatial(2,:)-spatial(1,:)))/length(spatial(1,:))

%plot spatial correlation
figure;plot(d/lambda,abs(spatial(1,:)),'b');
hold on;
plot(d/lambda,abs(spatial(2,:)),'r');
plot(d/lambda,abs(spatial_num),'g');
plot(d/lambda,abs(spatial_circle),'p');
axis([0 1 0 1]);
xlabel('Antenna Separation in wavelength');
ylabel('Spatial Correlation')
grid on;
%%
bitseq = randn(1,100);
bitseq(bitseq>0)=1;
bitseq(bitseq<0)=0;
QPSK_table = [1+1j  -1+1j -1-1j 1-1j]/sqrt(2); 
t_symbol = linspace(0,1,32);f = 1;
for i=1:length(bitseq)/2
   temp = bitseq(2*(i-1)+1)*2 +bitseq(2*(i-1)+2);
   QPSK_symbols(i) = QPSK_table(temp+1);
   sig((i-1)*32+1:i*32) = real(QPSK_symbols(i)).*cos(2*pi*f*t_symbol)-1j*imag(QPSK_symbols(i)).*sin(2*pi*f*t_symbol);
end

for i = 1:length(d)
    for j = 1:length(ideal_phi)
        h(i,1,j) = alpha*exp(1j*(beta)).*sqrt(PAS(j)/sum_PAS);
        h(i,2,j) = alpha*exp(1j*(beta + 2*pi*(d(i)/lambda)*sin(ideal_phi(j) - phi_a))).*sqrt(PAS(j)/sum_PAS);
        sig_rx(i,1,j,:) = conv(sig,h(i,1,j));
        sig_rx(i,2,j,:) = conv(sig,h(i,1,j));
    end
    
    
    for k = 1:length(real_phi)
        tau(k) = d(i)*sin(real_phi(k)-phi_a)/c;
        spatial(2,i) = spatial(2,i) + exp(-1i*2*pi*fc * (tau(k) + error_tk(k)) ).*real_PAS(k)/real_sum_PAS;
    end
end

% [a,b] = xcorr(h(2,:),h(1,:));
% dd = abs(imag(a));

% stat = sqrt(sum((abs(spatial(2,:)-spatial(1,:)).^2))/length(spatial(1,:)))
stat = sum(abs(spatial(2,:)-spatial(1,:)))/length(spatial(1,:))


for i = 1:length(d)
    spatial_num(i,:) = xcorr(reshape(sig_rx(i,1,:),1,length(sig_rx)),reshape(conj(sig_rx(i,2,:)),1,length(sig_rx)));
    imag_a(i,:) = abs(imag(spatial_num(i,:)));
    
    imag_a_max(i) = max(abs(imag(spatial_num(i,:))));
end

figure;plot(d/lambda,abs(spatial(1,:)),'b');
hold on;
plot(d/lambda,abs(spatial(2,:)),'r');
% axis([0 1 0 1]);
xlabel('Antenna Separation in wavelength');
ylabel('Spatial Correlation');
grid on;

%figure;surf(imag_a);
%title(scenario);
%figure;
plot(d/lambda,imag_a(:,180));
% figure;plot(imag_a_max);