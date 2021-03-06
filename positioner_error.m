%%??????????????????????????????????????????
clc;clear;close all;
v1=15:5:165;
v2=[13.25	18.15	23.05	27.95	32.95	38	42.95	47.9	52.85	58.05...
    63.1	68.05	73	77.85	82.8	87.75	92.6	97.6	102.6	107.5...
    112.4	117.35	122.35	127.2	132.2	136.95	142	146.95	151.9	156.9	162
    ];
%%??????????
v3=v2+(15-v2(1,1));%%??15�?????
error_v=v1-v3;
x=15:1:165;
error_v=interp1(v1,error_v,x,'linear');
error_v_compensation=cat(1,x,error_v);
erro_v_max=max(error_v);
p = polyfit(x,error_v,1);
x1 = 15:0.2:165;
y1 = polyval(p,x1);
figure;plot(x,error_v,'LineWidth',2);
hold on;plot(x1,y1,'LineWidth',2);
xlabel('Vertical angle(deg)');ylabel('Angle error');
legend('angle error','linear fitting');

%%??????????
h1=0:5:360;
h2=[0	4.8	9.85	14.9	19.85	24.8	29.9	34.8	39.8	44.75	49.7...
    54.7	59.7	64.6	69.55	74.6	79.5	84.65	89.6	94.55	99.45...
    104.45	109.45	114.4	119.35	124.45	129.4	134.45	139.45	144.35	149.5...
    154.5	159.6	164.5	169.45	174.4	179.6	184.65	189.5	194.35	199.45...
    204.5	209.6	214.55	219.6	224.45	229.6	234.65	239.6	244.75	249.75...
    254.75	259.7	264.7	269.8	274.8	279.8	284.85	289.9	294.95	299.9...
    304.8	310	314.9	320.05	324.9	330.05	334.95	339.9	344.95	350.05	354.95	359.95
    ];
y=0:1:360;
error_h=h1-h2;
error_h=interp1(h1,error_h,y,'linear');
error_h_compensation=cat(1,y,error_h);
erro_h_max=max(error_h);
figure;plot(y,error_h,'LineWidth',2);
xlabel('Horizontal angle(deg)');ylabel('Angle error');