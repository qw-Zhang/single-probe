%%
%generate Laplace distribution
x = linspace(-pi,pi,100);
miu = [0.7,-13.2,146.1,-30.5,-11.4,-1.1]*pi/180;b = 35*pi/180;
figure;hold on;
for i = 1:6
Lap(i,:) =(1/(2*b))*exp(-abs(x-miu(i))/b);
% integral(Lap,-pi,pi)
% plot(x,Lap(i,:));
end
pow = [-3,-4.3,-5.7,-7.3,-9,-11.4];
amp = 10.^(pow/10);
amp = amp./max(amp);
for i = 1:6
    LAP(i,:) = amp(i).*Lap(i,:);
    plot(x,LAP(i,:));
end
Lap_fin = sum(LAP); Lap_fin = Lap_fin/sum(Lap_fin);
%% 
fc = 2.535e8;
c = 3e8;
lambda = c/fc;
d = linspace(0,lambda,100);

phi_a = 0*pi/180;
scenario = 'macro';

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
%PAS continous->distrete    laplace->impulse
AOA = round([0.7,-13.2,146.1,-30.5,-11.4,-1.1]);
phi_a = 60*pi/180;
pow = [-3,-4.3,-5.7,-7.3,-9,-11.4];
amp = 10.^(pow/10);
A = amp./max(amp);
spatial = zeros(1,length(d));
PAS = zeros(1,360);
for i = 1:6
    PAS(AOA(i)+180) = A(i);
end
% PAS = PAS*pi/180;
sum_PAS = sum(PAS);
phi_p = linspace(-pi,pi,360);

for i = 1:100
    temp = 0;
    for j = 1:length(phi_p)
        temp = temp + exp(-1i*2*pi*(d(i)/lambda)*sin(phi_p(j) - phi_a)).*PAS(j)/sum_PAS;
    end
    spatial(i) = temp;
end

figure;plot(d/lambda,abs(spatial));
axis([0,1,0,1]);