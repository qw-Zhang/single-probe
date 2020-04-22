clc
clear
path_num=5;
% Ratio=0.7;
Attenuation=3;%db
near_point_num=9;
m=1;

delay1_time=0;%ms
delay1_point=floor( delay1_time/(18.75/768) );


delay2_time=2;%ms
delay2_point=floor( delay2_time/(18.75/768) );


delay3_time=5;%ms
delay3_point=floor( delay3_time/(18.75/768) );


delay4_time=10;%ms
delay4_point=floor( delay4_time/(18.75/768) );


delay5_time=15;%ms
delay5_point=floor( delay5_time/(18.75/768) );

if Attenuation==0;%db每径能量相等threhold=0.8
    amp1=1;
    amp2=1;
    amp3=1;
    amp4=1;
    amp5=1;
elseif Attenuation==3;%逐径衰减3db
    amp1=1;
    amp2=1/sqrt(2);
    amp3=1/sqrt(4);
    amp4=1/sqrt(8);
    amp5=1/sqrt(16);
elseif Attenuation==6;%逐径衰减6db
    amp1=1;
    amp2=1/2;
    amp3=1/4;
    amp4=1/8;
    amp5=1/16;
end 


sdandard_timeoffset=[delay1_point,delay2_point,delay3_point,delay4_point,delay5_point]


data_n=35;     % 数据段每个符号中的有用点数
data_N=128;     % 数据段每个符号的IFFT点数
cp_N=data_N/2;       %数据段循环前缀长度
x=8*4;
n=data_n*x;    % 同步段每个符号中的有用点数
N=data_N*x;     % 同步段每个符号的IFFT点数 4096

data_N=4096*1+3;
%=============================================
%=========get syn local data===================
% rand('seed',0);
%  for m=1:10
PN_n=exp(i*2*pi*rand(1,n));
temp=zeros(1,N);
temp([N-n/2+1:N,1:n/2])=PN_n;

time_data=ifft(temp);
enrg_localtime=sum(abs(time_data).^2);
local=(conj([fliplr(time_data),zeros(1,N)]));
local_freq=fft(local);
%=======同步头后补充一段随机数据
  %syn_data=[time_data,exp(i*2*pi*rand(1,data_N))/1];
 syn_data=[time_data,zeros(1,data_N)];
 tmp = rand(1,delay5_point);
path1_syn_data=amp1*syn_data;
path2_syn_data=amp2*[exp(i*2*pi*tmp(1:delay2_point))/1000,syn_data];
path3_syn_data=amp3*[exp(i*2*pi*tmp(1:delay3_point))/1000,syn_data];
path4_syn_data=amp4*[exp(i*2*pi*tmp(1:delay4_point))/1000,syn_data];
path5_syn_data=amp5*[exp(i*2*pi*tmp(1:delay5_point))/1000,syn_data];

L=length(path1_syn_data);

if path_num==1
rev_syn_data=path1_syn_data(1:L);
elseif path_num==2
rev_syn_data=path1_syn_data(1:L)+path2_syn_data(1:L);
elseif path_num==3
rev_syn_data=path1_syn_data(1:L)+path2_syn_data(1:L)+path3_syn_data(1:L);
elseif path_num==4
rev_syn_data=path1_syn_data(1:L)+path2_syn_data(1:L)+path3_syn_data(1:L)+path4_syn_data(1:L);
elseif path_num==5
rev_syn_data=path1_syn_data(1:L)+path2_syn_data(1:L)+path3_syn_data(1:L)+path4_syn_data(1:L)+path5_syn_data(1:L);

end
rev_syn_data_prezero=[zeros(1,N),rev_syn_data];


%========fft对应圆周卷积，以圆周卷积来实现线性卷积，同步所用的共轭相关算法可以用线性卷积实现，求得归一化相关值＝＝＝＝＝
u=1;
enrg_temp=0;
for k=2:N+1
enrg_temp=enrg_temp+abs(rev_syn_data_prezero(k)).^2;
end
for i=0:(length(rev_syn_data_prezero)/N)-2
    
     now_syn_data=rev_syn_data_prezero(i*N+1:(i+2)*N);
     syn_freq=fft( now_syn_data  );
     temp=syn_freq.*local_freq;
     temp_time=ifft(temp);
     corrl(i*N+1:(i+1)*N)=temp_time(N+1:2*N);
     amp_corrl(i*N+1:(i+1)*N)=abs(corrl(i*N+1:(i+1)*N)).^2;
     
     for k=1:N
         uni_amp_corrl(i*N+k)=amp_corrl(i*N+k)/enrg_temp/enrg_localtime;
         enrg_temp=enrg_temp-abs(rev_syn_data_prezero(i*N+k))^2;
         enrg_temp=enrg_temp+abs(rev_syn_data_prezero(i*N+k+N+1))^2;
     end
 end
 %=========门限百分比（最低峰值/最高峰值）threhold_est，用于搜索，以便最后确定一个最优门限=============================
 temp_timeoffset=sdandard_timeoffset+4096;
 max_uni_correl=uni_amp_corrl( floor (temp_timeoffset) ) ;
 threhold_est(m)=max_uni_correl(5)/max_uni_correl(1);
%  m=m+1;
%  end
% min(threhold_est)





 if Attenuation==0
   threhold_est=0.8;
 end
 %=========================================================================

 
 for i=0:(length(rev_syn_data_prezero)/N)-2
     for k=1:N
      if uni_amp_corrl(i*N+k)>=max(uni_amp_corrl)*threhold_est
         timeoffset(u)=i*N+k-N
         u=u+1;
      end
    end
 end
     

 path_num_est=u-1
%  delay1_time_est=timeoffset(1)*(18.75/768)
%  delay2_time_est=timeoffset(2)*(18.75/768)
%  delay3_time_est=timeoffset(3)*(18.75/768)
%  delay4_time_est=timeoffset(4)*(18.75/768)
% 
%   test_erng=sum(  abs(rev_syn_data_prezero(8194:8194+4095) ).^2 ) ;  
 figure(1)
 plot(uni_amp_corrl);
 
 
 for k=1:length(temp_timeoffset)
     near_ID=floor( temp_timeoffset(k)-near_point_num:temp_timeoffset(k)+near_point_num );
     near_max_uni_correl( (k-1)*(2*near_point_num+1)+1:k*(2*near_point_num+1) )=uni_amp_corrl (near_ID) ;
 end 
%  for k=1:length(temp_timeoffset)
%      near_ID=floor( temp_timeoffset(k)-near_point_num:temp_timeoffset(k)+near_point_num );
%     part_corrl=uni_amp_corrl (near_ID) 
%  end 
 figure(2)
 plot(near_max_uni_correl);
%===================由粗略估计的一组timeoffset进行精确估计＝＝＝＝＝＝＝＝＝＝＝＝＝＝＝ 
 precis1_timeoffset=timeoffset(find(timeoffset>-100));
 for k=1:length(precis1_timeoffset)
     
 end

%%
for i = 1:5
rr(1,i) = cursor_info(i).Position(2)./cursor_info(5).Position(2);
rr(2,i) = cursor_info1(i).Position(2)./cursor_info1(5).Position(2);
end
rr(1,:)
power(rr(1,:),2)
[amp1,amp2,amp3,amp4,amp5].^2