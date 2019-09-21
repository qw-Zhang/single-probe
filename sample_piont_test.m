f0 = 10;f1 = 2000;
t = linspace(0,9999,10000);
chirp = cos(2 * pi*(f0*t / 10000 + (f1 / 2)*power(t / 10000, 2))) + 1i*sin(2 * pi*(f0*t / 10000 + (f1 / 2)*power(t / 10000, 2)));

figure;plot(b_d,a_d);
rate = 100;
chirp_in = resample(chirp,rate,1);
chirp_ex_len = length(chirp_in)/rate;

piont = 11;

j = 1;
for i = 1:chirp_ex_len
    chirp_piont(j) = chirp_in(rate*(i-1)+piont);
    j = j + 1;
end

% for r = 1:rate
%     j = 1;
%     for i = 1:chirp_ex_len
%         chirp_ex(r,j) = chirp_in(rate*(i-1)+r);
%         j = j + 1;
%     end
%     [a,b] = xcorr(chirp_piont,chirp_ex(r,:));
%     [v_m(2,r),p_temp] = max(a);
%     v_m(1,r) = b(p_temp);
%     
% end
% figure;plot(v_m(2,:));
% [r_m(2),r_m(1)] = max(v_m(2,:));

chirp_piont_in = resample(chirp_piont,rate,1);
for r = 1:rate
    j = 1;
    for i = 1:chirp_ex_len
        chirp_piont_ex(r,j) = chirp_piont_in(rate*(i-1)+r);
        j = j + 1;
    end
    [a,b] = xcorr(chirp_piont_ex(r,:),chirp);
    [v_in_m(2,r),p_temp] = max(abs(a));
    v_in_m(1,r) = b(p_temp);
end
% figure;plot(v_in_m(2,:));
[r_m(4),r_m(3)] = max(v_in_m(2,:));
