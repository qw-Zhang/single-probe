function repair = IQrepair(data,tx)
    [a,b] = xcorr(real(data),tx);
    location = find_start(a);
    start = location - 11;
    p = start + 50;
    temp_real = real(data(p));
    temp_imag = imag(data(p));
    temp1 = abs(temp_real - temp_imag);
    temp2 = abs(temp_real + temp_imag);
    if temp1 > temp2
        repair = 1;
    else
        repair = 0;
    end
end