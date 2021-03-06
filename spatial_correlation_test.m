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
grid on;

%%
%PAS continous->distrete
clear;
fc = 2.535e9;
c = 3e8;
lambda = c/fc;
d = linspace(0,lambda,100);
ideal_phi = linspace(-pi,pi,360);
PAS = zeros(1,length(ideal_phi));

real_phi = linspace(-pi,pi,36);
real_PAS = zeros(1,length(real_phi));

rate = length(ideal_phi)/length(real_phi);

phi_a = 0*pi/180;
scenario = 'micro';

switch scenario
    case 'test'
        AOA = -150*pi/180;AS = 35*pi/180;
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

% ideal scenario
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
% error_tk = 1e-12*randn(1,length(real_phi)); %the error of synthesis STSs
error_tk = zeros(1,length(real_phi));
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
figure;plot(ideal_phi*180/pi,PAS)
figure;plot(d/lambda,abs(spatial(1,:)),'b');
hold on;
plot(d/lambda,abs(spatial(2,:)),'r');
axis([0 1 0 1]);
xlabel('Antenna Separation in wavelength');
ylabel('Spatial Correlation')
grid on;

%%
%test with h1 h2 without the euqalation of saptial correlation
%using three method generate spatial correlation:
%   1. traditional formula method
%   2. numberical calculation method
%   3. simulation two antennas with circle method
%   4. add signal transfer based method_3
clear;
fc = 2.535e9;
c = 3e8;
lambda = c/fc;
d = linspace(0,lambda,100);
new_d = linspace(0,lambda/2,100);
ideal_phi = linspace(-pi,pi,360);

real_phi = linspace(-pi,pi,18);
real_PAS = zeros(1,length(real_phi));

rate = length(ideal_phi)/length(real_phi);

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
        error('wrong scenario.');
end

% %ideal scenario
% amp = 10.^(pow/10);
% A = amp./max(amp);
% for i = 1:length(ideal_phi)
%     for j = 1:length(AOA)
%         PAS(i) = PAS(i) + A(j)*(1/(2*AS))*(exp(-abs(ideal_phi(i)-AOA(j))/AS));
%     end
% end

%ideal scenario
amp = 10.^(pow/10);
A = amp./max(amp);
for i = 1:length(ideal_phi)
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

%real probe scenario
for i = 1:length(real_phi)
    real_PAS(i) = sum(PAS(rate*(i-1)+1:rate*i))/rate;
    real_phi(i) = sum(ideal_phi(rate*(i-1)+1:rate*i))/rate;
end


%Rx one antenna pattern
ant = design(dipole ,fc);
%generate the pattern gain of the antenna.
%parameter:antenna,frequency,al,ez.
%G(dBi)=10lgGi G(dBd)=10lgG dBi=dBd+2.15
% [D,al,ez] = pattern(ant,fc,0:1:360,-90:1:90);
[P_rx] = patternAzimuth(ant,fc,'Azimuth',-180:1:180);

%Define Rx Array
arrayObject = linearArray;
arrayObject.NumElements = 2;
%Define Array Elements
Element1 = dipole;
arrayObject.Element = [Element1 ];
% Design array at frequency 2535000000Hz
arrayObject = design(arrayObject,2535000000,Element1);
%arrayObject Properties Changed
arrayObject.ElementSpacing = 0.05;
% Azimuth for linearArray
plotFrequency = 2535000000; azRange = 0:1:360; Termination = 50;
figure;
[P_az] = pattern(arrayObject, plotFrequency,azRange,0,'Termination',Termination);

spatial = zeros(2,length(d));
beta = -pi + (pi--pi).*rand(1,length(ideal_phi));
% beta = 0;
alpha = 1;

%using for normalized PAS
sum_PAS = sum(PAS);
real_sum_PAS = sum(real_PAS);

%random error
% !!!   1/100 nanosecond...
% error_tk = 1e-11*randn(1,length(real_phi));
error_tk = zeros(1,length(ideal_phi));


%using numberical method calculate the spatial correlation
spatial_num = zeros(length(d),1);

%simulate radius of circle
r = 1;

t = linspace(0,1e-5,100);
% doppler frequency fd = fc*v*cos(theta)/c
v = 10; %m/s
theta = 30*pi/180;
fd = fc*v*cos(theta)/c;
sig = exp(1i*2*pi*(fc + fd)*t);
rx_1 = repmat(1+1j,length(d),length(ideal_phi),length(sig));
rx_2 = repmat(1+1j,length(d),length(ideal_phi),length(sig));
spatial_circle_sig = zeros(1,length(d));
spatial_circle = zeros(1,length(d));
[d_1, d_2, delta_d, h_sig_1, h_sig_2] = deal(zeros(length(d), length(ideal_phi)));
for i = 1:length(d)
    %simulate spatial correlation using two antennas with circle
    for j = 1:length(ideal_phi)
        phi_1 = (ideal_phi(j) - phi_a);
        phi_2 = pi - phi_1;
        d_1(i,j) = sqrt( r^2 + new_d(i)^2 - 2*new_d(i)*r*cos( phi_1 + pi/2 ) );
        d_2(i,j) = sqrt( r^2 + new_d(i)^2 - 2*new_d(i)*r*cos( phi_2 - pi/2 ) );
        delta_d(i,j) = d_2(i,j) - d_1(i,j);
        h_sig_1(i,j) = alpha*exp(1j*(beta(j))).*sqrt(PAS(j)/sum_PAS);
        h_sig_2(i,j) = alpha*exp(1j*(beta(j) + 2*pi*delta_d(i,j)/lambda)).*sqrt(PAS(j)/sum_PAS);
        rx_1(i,j,:) = conv(h_sig_1(i,j),sig);
        rx_2(i,j,:) = conv(h_sig_2(i,j),sig);
    end
    
    %     Corr(:,:,i) = abs(corrcoef(h_new_1(i,:),h_new_2(i,:)));
    %     spatial_circle_sig(i) = sum(r1(i,:).*conj(r2(i,:)))/sum(r1(1,:).*conj(r2(1,:)));
    
    num = randi([1,length(d)],1);
    %randomly choose one point in signal sequence. "100" is chosen randomly
    spatial_circle_sig(i) = sum(rx_1(i,:,num).*conj(rx_2(i,:,num)));
    spatial_circle(i) = sum(h_sig_1(i,:).*conj(h_sig_2(i,:)));
    
    %simulate spatial correlation using two antennas with circle, 
    %real phi and PAS
    
    for j = 1:length(real_phi)
        phi_real_1 = (real_phi(j) - phi_a);
        phi_real_2 = pi - phi_real_1;
        d_real(1,i,j) = sqrt( r^2 + new_d(i)^2 - 2*new_d(i)*r*cos( phi_real_1 + pi/2 ) );
        d_real(2,i,j) = sqrt( r^2 + new_d(i)^2 - 2*new_d(i)*r*cos( phi_real_2 - pi/2 ) );
        delta_real_d(i,j) = d_real(2,i,j) - d_real(1,i,j);
        h_sig_real_1(i,j) = alpha*exp(1j*(beta(j))).*sqrt(real_PAS(j)/real_sum_PAS);
        h_sig_real_2(i,j) = alpha*exp(1j*(beta(j) + 2*pi*delta_real_d(i,j)/lambda)).*sqrt(real_PAS(j)/real_sum_PAS);
        rx_real_1(i,j,:) = conv(h_sig_real_1(i,j),sig);
        rx_real_2(i,j,:) = conv(h_sig_real_2(i,j),sig);
    end
    
    %     Corr(:,:,i) = abs(corrcoef(h_new_1(i,:),h_new_2(i,:)));
    %     spatial_circle_sig(i) = sum(r1(i,:).*conj(r2(i,:)))/sum(r1(1,:).*conj(r2(1,:)));
    
    num = randi([1,length(d)],1);
    %randomly choose one point in signal sequence. "100" is chosen randomly
    spatial_circle_real_sig(i) = sum(rx_real_1(i,:,num).*conj(rx_real_2(i,:,num)));
    spatial_circle_real(i) = sum(h_sig_real_1(i,:).*conj(h_sig_real_2(i,:)));
    
    beta1 = 0;
    %using numberical method calculate the spatial correlation
    for j = 1:length(ideal_phi)
        h(i,1,j) = alpha*exp(1j*(beta1)).*sqrt(PAS(j)/sum_PAS);
        h(i,2,j) = alpha*exp(1j*(beta1 + 2*pi*(2*new_d(i)/lambda) * cos((pi/2-ideal_phi(j)) + phi_a) ) ).*sqrt(PAS(j)/sum_PAS);
    end
    h1_cal = reshape(h(i,1,:),1,length(PAS));
    h2_cal = reshape(conj(h(i,2,:)),1,length(PAS));
    spatial_num(i,:) = sum(h1_cal.*h2_cal);
    
    Corr(i,:,:) = abs(corrcoef(squeeze(h(i,1,:)),squeeze(h(i,2,:))));
    
    %using traditional method(ideally equation) calculate
    for k = 1:length(ideal_phi)
        tau(k) = d(i)*sin(ideal_phi(k)-phi_a)/c;
        spatial(2,i) = spatial(2,i) + exp(-1i*2*pi*fc * (tau(k) + error_tk(k)) ).*PAS(k)/sum_PAS;
    end
end

Correaltion = squeeze(Corr(:,1,2));

%calculate the statistical error
stat = sum(abs(spatial_circle_real_sig-spatial_num'))/length(spatial_num)

%plot spatial correlation
figure;
%plot(d/lambda,abs(spatial(1,:)),'b');
hold on;
% plot(d/lambda,abs(spatial(2,:)),'r');
plot(d/lambda,abs(spatial_num),'green');
plot(d/lambda,abs(spatial_circle_real_sig),'red');
plot(d/lambda,abs(spatial_circle_sig),'blue');
% plot(d/lambda,Correaltion,'p');
% axis([0 1 0 1]);
xlabel('Antenna Separation in wavelength');
ylabel('Spatial Correlation');
grid on;
% legend('theo equation','theo num','sim circle h','sim circle sig');
legend('theo num','sim circle real sig','sim circle sig');