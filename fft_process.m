function sig_fft = fft_process(sig)
fft_temp = fft(sig);
[val,pos] = max(abs(fft_temp));
sig_fft = fft_temp(pos);
