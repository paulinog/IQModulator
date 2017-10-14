figure('Position', [26 340 1520 420]),
paramFilt.SampleRate = SampleRate;
paramFilt.BW = 1e3;
paramFilt.freq_central = f_I;
paramFilt.order = 2;
paramFilt.gain = 1;
paramFilt.plot_flag = false;
FilteredToneI = custom_filter(PD_Signal, paramFilt);
Spectrum_FilteredToneI = abs(fftshift(fft(FilteredToneI)));
subplot(251),plot(time,FilteredToneI)
title('Filtered Pilot Tone I'), xlabel('Samples')
legend(['f_I=' num2str(f_I) ' Hz'])
subplot(256),plot(frequency,10*log10(Spectrum_FilteredToneI),'k')
title('Spectrum of the Filtered Pilot tone I'), axis([f_I-1e3 f_I+1e3 -60 inf])
xlabel('Frequency'), ylabel('Power Spectral Density'), ylim([-20 25])

paramFilt.SampleRate = SampleRate;
paramFilt.BW = 1e3;
paramFilt.freq_central = f_Q;
paramFilt.order = 2;
paramFilt.gain = 1;
paramFilt.plot_flag = false;
FilteredToneQ = custom_filter(PD_Signal, paramFilt);
Spectrum_FilteredToneQ = abs(fftshift(fft(FilteredToneQ)));
subplot(252),plot(time,FilteredToneQ)
title('Filtered Pilot Tone Q'), xlabel('Samples')
legend(['f_Q=' num2str(f_Q) ' Hz'])
subplot(257),plot(frequency,10*log10(Spectrum_FilteredToneQ),'k')
title('Spectrum of the Filtered Pilot tone Q'), axis([f_Q-1e3 f_Q+1e3 -60 inf])
xlabel('Frequency'), ylabel('Power Spectral Density'), ylim([-20 25])


paramFilt.SampleRate = SampleRate;
paramFilt.BW = 1e3;
paramFilt.freq_central = 2*f_I;
paramFilt.order = 2;
paramFilt.gain = 1;
paramFilt.plot_flag = false;
FilteredTone2I = custom_filter(PD_Signal, paramFilt);
Spectrum_FilteredTone2I = abs(fftshift(fft(FilteredTone2I)));
subplot(253),plot(time,FilteredTone2I)
title('2nd Harmonic of I'), xlabel('Samples')
legend(['2*f_I=' num2str(2*f_I) ' Hz'])
subplot(258),plot(frequency,10*log10(Spectrum_FilteredTone2I),'k')
title('Spectrum of 2nd Harmonic of I'), axis([2*f_I-1e3 2*f_I+1e3 -60 Inf])
xlabel('Frequency'), ylabel('Power Spectral Density'), ylim([-20 25])

paramFilt.SampleRate = SampleRate;
paramFilt.BW = 1e3;
paramFilt.freq_central = 2*f_Q;
paramFilt.order = 2;
paramFilt.gain = 1;
paramFilt.plot_flag = false;
FilteredTone2Q = custom_filter(PD_Signal, paramFilt);
Spectrum_FilteredTone2Q = abs(fftshift(fft(FilteredTone2Q)));
subplot(254),plot(time,FilteredTone2Q)
title('2nd Harmonic of Q'), xlabel('Samples')
legend(['2*f_Q=' num2str(2*f_Q) ' Hz'])
subplot(2,5,9),plot(frequency,10*log10(Spectrum_FilteredTone2Q),'k')
title('Spectrum of 2nd Harmonic of Q'), axis([2*f_Q-1e3 2*f_Q+1e3 -60 Inf])
xlabel('Frequency'), ylabel('Power Spectral Density'), ylim([-20 25])

paramFilt.SampleRate = SampleRate;
paramFilt.BW = 1e3;
paramFilt.freq_central = f_Q-f_I;
paramFilt.order = 2;
paramFilt.gain = 1;
paramFilt.plot_flag = false;
FilteredBeating = custom_filter(PD_Signal, paramFilt);
Spectrum_FilteredBeating = abs(fftshift(fft(FilteredBeating)));
subplot(255),plot(time,FilteredBeating)
title('Filtered Lower Beating'), xlabel('Samples')
legend(['f_Q-f_I=' num2str(f_Q-f_I) ' Hz'])
subplot(2,5,10),plot(frequency,10*log10(Spectrum_FilteredBeating),'k')
title('Spectrum of the Lower Beating'), axis([(f_Q-f_I-1e3) (f_Q-f_I+1e3) -60 Inf])
xlabel('Frequency'), ylabel('Power Spectral Density')