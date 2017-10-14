function filtered_signal_time = custom_filter(signal,param)
    [row,col] = size(signal);
    frequency = linspace(-0.5, 0.5,col)*param.SampleRate;
    filtered_signal_time = zeros(row,col);
    
    GAIN = param.gain;
    
    filt = exp(-1*log(sqrt(2))*((frequency - param.freq_central)/(param.BW/2)).^(2*param.order));
    filt = filt + fliplr(filt);
    filt = filt/max(filt);
    
    for i=1:row
        filtered_signal_freq = fftshift(fft(signal(i,:))).*filt;
        filtered_signal_time(i,:) =  GAIN*real(ifft(ifftshift(filtered_signal_freq)));
     end
            
    if param.plot_flag
        figure('Position',[318 342 900 280]),
        subplot(131),plot(frequency,abs(fftshift(fft(signal))))
        title('FFT of the Signal'), xlabel('Frequency')
        subplot(132),plot(frequency,filt)
         title('Filter Response'), xlabel('Frequency')
        subplot(133),plot(frequency,abs(filtered_signal_freq))
         title('FFT of the Filtered Signal'), xlabel('Frequency')
    end
end

