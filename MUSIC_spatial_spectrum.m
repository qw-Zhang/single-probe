% clear; close all;
%%%%%%%% MUSIC for Uniform Linear Array%%%%%%%%
derad = pi/180;
N = 8;               % 阵元个数        
M = 3;               % 信源数目
theta = [-30 0 60];  % 待估计角度
snr = 10;            % 信噪比
K = 512;             % 快拍数
 
dd = 0.5;            % 阵元间距 
d=0:dd:(N-1)*dd;
A=exp(-1i*2*pi*d.'*sin(theta*derad));
% 构建信号模型
S=randn(M,K); X=A*S;
X1=awgn(X,snr,'measured');
% 计算协方差矩阵
Rxx=X1*X1'/K;
% 特征值分解
[EV,D]=eig(Rxx);
EVA=diag(D)';
[EVA,I]=sort(EVA);
EV=fliplr(EV(:,I));
En=EV(:,M+1:N);    % 信号子空间
% 遍历每个角度，计算空间谱
for iang = 1:361
    angle(iang)=(iang-181)/2;
    phim=derad*angle(iang);
    a=exp(-1i*2*pi*d*sin(phim)).';    
    Pmusic(iang)=1/(a'*En*En'*a);
end
Pmusic=abs(Pmusic);
Pmusic=10*log10(Pmusic);
h=plot(angle,Pmusic);
set(h,'Linewidth',0.5);
xlabel('入射角/(degree)');
ylabel('空间谱/(dB)');
set(gca, 'XTick',[-90:30:90]);
grid on;