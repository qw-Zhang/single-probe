% clear; close all;
%%%%%%%% MUSIC for Uniform Linear Array%%%%%%%%
derad = pi/180;
N = 8;               % ��Ԫ����        
M = 3;               % ��Դ��Ŀ
theta = [-30 0 60];  % �����ƽǶ�
snr = 10;            % �����
K = 512;             % ������
 
dd = 0.5;            % ��Ԫ��� 
d=0:dd:(N-1)*dd;
A=exp(-1i*2*pi*d.'*sin(theta*derad));
% �����ź�ģ��
S=randn(M,K); X=A*S;
X1=awgn(X,snr,'measured');
% ����Э�������
Rxx=X1*X1'/K;
% ����ֵ�ֽ�
[EV,D]=eig(Rxx);
EVA=diag(D)';
[EVA,I]=sort(EVA);
EV=fliplr(EV(:,I));
En=EV(:,M+1:N);    % �ź��ӿռ�
% ����ÿ���Ƕȣ�����ռ���
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
xlabel('�����/(degree)');
ylabel('�ռ���/(dB)');
set(gca, 'XTick',[-90:30:90]);
grid on;