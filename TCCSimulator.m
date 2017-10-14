%% IQ Modulator Bias Control - Simulator
% Guilherme Paulino
%
% Trabalho de Conclusão de Curso (Final Design)
%
% University of Campinas, UNICAMP. 2017
%
% Version 0.1 (local)

%% Initialization
close all; clear all;

%% IQ Modulator Parameters
Vpi.I = 3;              % Half-wave voltage, Typ: 3.75 Vpp
Vpi.Q = 3;
Vpi.P = 3;
IL = 0;                 % Insertion Loss,  Typ: 8.5 dB
ER_I = 10^( 30 /10);    % Extinction Ratio, Typ: 30 dB
ER_Q = 10^( 30 /10);

% Sweep ----------> move to the right place
SweepPts = 300;
SweepRangeMax = 3*Vpi.I;

% Signal parameters
SymbolRate = 1e6;                           % 1 MHz
SamplePerSymbol = 4;
SampleRate = SymbolRate*SamplePerSymbol;    % 4 MSa/s
FrequencyResolution = 1e2;                  % 100 Hz
L = SampleRate/FrequencyResolution;         % vector length

% Time/Frequency domain vectors
time = 0 : 1/(SampleRate-1) : (1/FrequencyResolution);
frequency = linspace(-0.5,0.5,L)*SampleRate;

%% PRBS Signal Generation
% * upsample(X,N) upsamples input signal X by inserting
%     N-1 zeros between input samples.  X may be a vector
%     or a signal matrix (one signal per column).
% * Signum function: For each element of X, sign(X) returns 1
%     if the element is greater than zero, 0 if it equals zero
%     and -1 if it is less than zero.
%     For the nonzero elements of complex X, sign(X) = X ./ ABS(X).
% * randn Normally distributed pseudorandom numbers.
%     R = randn(N) returns an N-by-N matrix containing pseudorandom
%     values drawn from the standard normal distribution.
%     randn(M,N) or randn([M,N]) returns an M-by-N matrix. randn(M,N,P,...)
%     or randn([M,N,P,...]) returns an M-by-N-by-P-by-... array.
%     randn returns a scalar.  randn(SIZE(A)) returns an array the
%     same size as A.
%
Symbol_L = L/SamplePerSymbol; % Symbol length
% Signal Generation
% Pseudorandom binary sequence (PRBS), in format Non-return-to-zero (NRZ),
% After up converter (Sample Rate = 4 * Symbol Rate)
PRBS_I = upsample(sign(randn(1, Symbol_L)), SamplePerSymbol);
PRBS_Q = upsample(sign(randn(1, Symbol_L)), SamplePerSymbol);
Signal_I = zeros(1,L);
Signal_Q = zeros(1,L);
for i=0:(SamplePerSymbol-1)
    Signal_I = Signal_I + circshift(PRBS_I,[0 i]);
    Signal_Q = Signal_Q + circshift(PRBS_Q,[0 i]);
end

figure('Position',[618 412 300 280]),
plot(Signal_I, Signal_Q, 'o-' ), xlabel('Signal I'), ylabel('Signal Q')
title('QPSK Generated Signal (Normalized)')

%% BW limitation
% Transfer function for the filter response:
%
% $$ TF_{BW}(f) = \exp(-\log(\sqrt{2})(\frac{f}{0.75*2})^2) $$
%
if false
    figure,
    subplot(211),plot(time,Signal_I)
    title('Signal I'), ylabel('Samples'), axis([0 1e-4 -1.1 1.1])
    subplot(212),plot(time,Signal_Q)
    title('Signal Q'), ylabel('Samples'), axis([0 1e-4 -1.1 1.1])
    
    paramRF.SampleRate = SampleRate/4;
    paramRF.BW = 0.75*SymbolRate;
    paramRF.freq_central = 0;
    paramRF.order = 1;
    paramRF.gain = 1;
    paramRF.plot_flag = true;
    Signal_I = custom_filter(Signal_I,paramRF);
    title('Filtered Signal I'), ylabel('Frequency')
    Signal_Q = custom_filter(Signal_Q,paramRF);
    title('Filtered Signal Q'), ylabel('Frequency')
    
    figure,
    subplot(211),plot(time,Signal_I)
    title('Filtered Signal I'), ylabel('Samples'), axis([0 1e-4 -1.1 1.1])
    
    subplot(212),plot(time,Signal_Q)
    title('Filtered Signal Q'), ylabel('Samples'), axis([0 1e-4 -1.1 1.1])
    
    figure,
    plot(Signal_I, Signal_Q, 'o-' ), xlabel('Signal I'), ylabel('Signal Q')
    title('QPSK Generated Signal (Normalized)')
end

%% Pilot Tone Generation
%
% $$V_{pilot\ tone} = \frac{V_{pp}}{2}V_{\pi}\cos(2 \pi f t + \theta)$$
%
% * fftshift Shift zero-frequency component to center of spectrum.

f_I = 2e3; % 2 kHz
f_Q = 3e3; % 3 kHz
f_P = 0e3; % no tone applied
theta_I = 0; % initial tone phase
theta_Q = 0;
theta_P = 0;
Vpp_I = 0.1; % tone amplitude (peak-to-peak)
Vpp_Q = 0.1;
Vpp_P = 0.0;
tone.I = (Vpp_I/2)*Vpi.I*cos(2*pi*f_I*time + theta_I); % signal
tone.Q = (Vpp_Q/2)*Vpi.Q*cos(2*pi*f_Q*time + theta_Q);
tone.P = (Vpp_P/2)*Vpi.P*cos(2*pi*f_P*time + theta_P);

figure,
subplot(211), plot(time,tone.I)
title('Pilot tone I')
subplot(212), plot(time,tone.Q)
title('Pilot tone Q')
figure,
subplot(211), plot(frequency,abs(fftshift(fft(tone.I))),'r')
title('Spectrum of Pilot tone I'), axis([0 1e4 -inf inf])
legend(['f_I = ' num2str(f_I) ' Hz'])
subplot(212), plot(frequency,abs(fftshift(fft(tone.Q))),'r')
title('Spectrum of Pilot tone Q'), axis([0 1e4 -inf inf])
legend(['f_Q = ' num2str(f_Q) ' Hz'])


%% Generate RF Signal
Noise.I = 10^( -10 /10) * randn(1, L); % Noise pattern -10 dB
Noise.Q = 10^( -10 /10) * randn(1, L);

RF.I = 1 * Vpi.I * Signal_I + Noise.I; % RF Signal
RF.Q = 1 * Vpi.Q * Signal_Q + Noise.Q;

figure('Position',[618 412 300 280]),
plot(RF.I, RF.Q, '.' ), xlabel('Samples')
title('RF Generated Signal')

%% Mach-Zehnder Modulator: Model 1
% Single ended RF; dual drive bias; same Vpi for RF and bias;
% ER infinite
%
% $$ E_{OUT} = E_{IN}\cos[(RF + tone + bias)\frac{\pi}{2 V_{\pi}}] $$

mzm1 = @(E_in, RF, tone, bias, Vpi) E_in * (cos((RF + tone + bias) * pi/(2*Vpi)));

%% Model 2
% Single ended RF; dual drive bias; same Vpi for RF and bias;
% ER finite
%
% $$ E_{OUT} = E_{IN}\{\cos[(RF + tone +
%   bias)\frac{\pi}{2V_{\pi}}]+i \sqrt{ \frac{1}{ER} } \sin[(RF + tone +
%   bias)\frac{\pi}{2V_{\pi}}]\} $$

mzm2 = @(E_in, RF, tone, bias, Vpi, ER) E_in * (cos((RF + tone + bias) * pi/(2*Vpi)) + ...
    1i*sqrt(1/ER)*sin((RF + tone + bias)*pi/(2*Vpi)));

%% Quadrature Function
% $$ E_{out,Q} = E_{out} \exp[(bias.P + tone.P)\frac{i \pi}{V_{\pi,P}}] $$

quad = @(E_OUT, bias) exp(1i*((bias + tone.P)*pi)./Vpi.P) .* E_OUT;

%% QPSK Modulation
E_IN = 1; % Magnitude of the input field

% Bias Voltage Point
bias.I = 1.0   * Vpi.I; % Optimized for QPSK Modulation
bias.Q = 1.0   * Vpi.Q;
bias.P = 0.5 * Vpi.P;

E_OUT_I = mzm1(E_IN, RF.I, tone.I, bias.I, Vpi.I);
E_OUT_Q = mzm1(E_IN, RF.Q, tone.Q, bias.Q, Vpi.Q);
E_OUT_IQ = E_OUT_I + quad(E_OUT_Q, bias.P);

figure('Position',[618 412 300 280]),
plot(E_OUT_IQ, '.'), xlabel('Real part'), ylabel('Imaginary part')
title('QPSK Optical Signal - Model 1'), axis([-1.5 1.5 -1.5 1.5])

%% Extinction Ratio
figure('Position',[488 342 560 280]),
ER_I = 10^( 10 /10); % Low Extinction Ratio
E_OUT_I = mzm2(E_IN, RF.I, tone.I, bias.I, Vpi.I, ER_I);
subplot(121)
plot(E_OUT_I, 'o'), ylim([-0.2 0.2])
title('Lower Extinction Ratio'), ylabel('Polarization phase'), xlabel('E field magnitude')
legend('ER = 10 dB')

ER_I = 10^( 30 /10); % Extinction Ratio, Typ: 30 dB, higher is better
E_OUT_I = mzm2(E_IN, RF.I, tone.I, bias.I, Vpi.I, ER_I);
subplot(122)
plot(E_OUT_I, 'o'), ylim([-0.2 0.2])
title('Higher Extinction Ratio'), ylabel('Polarization phase'), xlabel('E field magnitude')
legend('ER = 30 dB')

%% QPSK Modulation, Model 2
ER_I = 10^( 20 /10); % Extinction Ratio, Typ: 30 dB, higher is better
ER_Q = 10^( 20 /10); %                   Min: 20 dB

E_OUT_I = mzm2(E_IN, RF.I, tone.I, bias.I, Vpi.I, ER_I);
E_OUT_Q = mzm2(E_IN, RF.Q, tone.Q, bias.Q, Vpi.Q, ER_Q);

bias.P = 0.5 * Vpi.P; % Optimum phase point
E_OUT_IQ = E_OUT_I + quad(E_OUT_Q, bias.P);

figure('Position',[618 412 300 280]),
plot(E_OUT_IQ, '.'), xlabel('Real part'), ylabel('Imaginary part')
title('QPSK Optical Signal - Model 2')
axis([-1.5 1.5 -1.5 1.5])

%%

% [param,Signal_I,Signal_Q] = pushbutton_GenerateRF_Callback(1,1,handles);
%
% param.Vpi_I = str2double(get(handles.edit_Vpi_I,'String'));
% param.Vpi_Q = str2double(get(handles.edit_Vpi_Q,'String'));
% param.Vpi_P = str2double(get(handles.edit_Vpi_P,'String'));
% param.ER_I = str2double(get(handles.edit_ExtRatio_I,'String'));
% param.ER_Q = str2double(get(handles.edit_ExtRatio_Q,'String'));
% param.ER_P = str2double(get(handles.edit_ExtRatio_P,'String'));
% param.IL = str2double(get(handles.edit_IL,'String'));
% contents = get(handles.popupmenu_Model,'String');
% param.Model = contents{get(handles.popupmenu_Model,'Value')};
%
% E_IN = 1;
%
% RF.I = Signal_I;
% RF.Q = Signal_Q;
%
% param.f_I = str2double(get(handles.edit_tone_freq_I,'String'));
% param.f_Q = str2double(get(handles.edit_tone_freq_Q,'String'));
% param.f_P = str2double(get(handles.edit_tone_freq_P,'String'));
% theta_I = str2double(get(handles.edit_tone_phase_I,'String'));
% theta_Q = str2double(get(handles.edit_tone_phase_Q,'String'));
% theta_P = str2double(get(handles.edit_tone_phase_P,'String'));
% Vpp_I = str2double(get(handles.edit_tone_Vpp_I,'String'));
% Vpp_Q = str2double(get(handles.edit_tone_Vpp_Q,'String'));
% Vpp_P = str2double(get(handles.edit_tone_Vpp_P,'String'));
% tone.I = (Vpp_I/2)*param.Vpi_I*cos(2*pi*param.f_I*param.time + theta_I);
% tone.Q = (Vpp_Q/2)*param.Vpi_Q*cos(2*pi*param.f_Q*param.time + theta_Q);
% tone.P = (Vpp_P/2)*param.Vpi_P*cos(2*pi*param.f_P*param.time + theta_P);
%
% bias.I = str2double(get(handles.edit_bias_I,'String'))*param.Vpi_I;
% bias.Q = str2double(get(handles.edit_bias_Q,'String'))*param.Vpi_Q;
% bias.P = str2double(get(handles.edit_bias_P,'String'))*param.Vpi_P;
%
% [E_OUT_IQ,E_OUT_I,E_OUT_Q] = modulate(E_IN,RF,bias,tone,param);
%
% plot(handles.axes_complex_diagram,real(E_OUT_I),imag(E_OUT_I),'*b')
% plot(handles.axes_complex_diagram,real(E_OUT_Q),imag(E_OUT_Q),'*r')
% plot(handles.axes_complex_diagram,real(E_OUT_IQ),imag(E_OUT_IQ),'*m'),xlim(handles.axes_complex_diagram,[-1 1]),ylim(handles.axes_complex_diagram,[-1 1])
%
% PD.R = str2double(get(handles.edit_PD_Responsivity,'String'))
% PD.Var = str2double(get(handles.edit_PD_AWGN,'String'))
% Signal_PD = monitor_PD(E_OUT_IQ,PD);
% Spectrum_PD = abs(fftshift(fft(Signal_PD)));
%
% if get(handles.checkbox_filter,'Value') == 1
%     get(handles.checkbox_filter,'Value')
%     paramFilt.SampleRate = param.SampleRate;
%     paramFilt.BW = str2double(get(handles.edit_BW,'String'));
%     paramFilt.freq_central = str2double(get(handles.edit_central_freq,'String'));
%     paramFilt.order = str2double(get(handles.edit_order,'String'));
%     paramFilt.gain = 1;
%     paramFilt.plot_flag = false;
%
%     Filtered = custom_filter(Signal_PD,paramFilt);
%     Spectrum_Filtered = abs(fftshift(fft(Filtered)));
%     plot(handles.axes_filtered_spectrum,param.frequency,20*log10(Spectrum_Filtered),'k'),xlim(handles.axes_filtered_spectrum,[-5e6 5e6]),ylim(handles.axes_filtered_spectrum,[-100 100])
% else
%     get(handles.checkbox_filter,'Value')
%     plot(handles.axes_filtered_spectrum,param.frequency,20*log10(Spectrum_PD),'r'),xlim(handles.axes_filtered_spectrum,[0 10e3]),ylim(handles.axes_filtered_spectrum,[-100 100])
%
% end

%% Photodetector - Average optical power
PD.Var = 10^(-90/10); % -90 dB of Noise floor
PD.R = 0.95;     % Responsivity
PD_Power = PD.R*abs(E_OUT_IQ).^2;
PD_Noise_floor = PD.Var*randn(size(E_OUT_IQ));

figure,plot(frequency, 10*log10(abs(fftshift(fft(PD_Noise_floor)))))
title('PD Noise Floor'), ylabel('Noise Power Spectral Density'), xlabel('Frequency')

% PD_Signal = E_OUT_IQ;           % TEST 1: Entry Signal
PD_Signal = PD_Power;           % TEST 2: Average Power
% PD_Signal = PD_Noise_floor;     % TEST 3: Noise
% PD_Signal = PD_Power + PD_Noise_floor; % Photodetector Signal
PD_Spectrum = abs(fftshift(fft(PD_Signal)));

figure, plot(time, 10*log10(PD_Signal)), xlabel('Samples')
title('Photodectector Received Signal')

figure,plot(frequency, 10*log10(PD_Spectrum),'r')
xlabel('Frequency'), ylabel('Power Spectral Density')
title('Spectrum of Photodectector Received Signal'), axis([0 1e4 -10 inf])

%% Heaater Transfer Function
paramFilt.SampleRate = SampleRate;
paramFilt.BW = 1e4;
paramFilt.freq_central = 0;
paramFilt.order = 1;
paramFilt.gain = 1;
paramFilt.plot_flag = false;

PD_Signal = custom_filter(PD_Signal, paramFilt);
Spectrum_PD_Signal = abs(fftshift(fft(PD_Signal)));

figure,plot(time,PD_Signal)
figure,plot(frequency,10*log10(Spectrum_PD_Signal),'k')
title('Heater Transfer Function'), axis([0 5e4 -60 Inf])

%%
paramFilt.SampleRate = SampleRate;
paramFilt.BW = 1e3;
paramFilt.freq_central = f_I;
paramFilt.order = 2;
paramFilt.gain = 1;
paramFilt.plot_flag = false;

FilteredToneI = custom_filter(PD_Signal, paramFilt);
Spectrum_FilteredToneI = abs(fftshift(fft(FilteredToneI)));

figure,plot(time,FilteredToneI)
figure,plot(frequency,10*log10(Spectrum_FilteredToneI),'k')
title('Spectrum of the Filtered Pilot tone I'), axis([0 1e4 -60 Inf])

%%
paramFilt.SampleRate = SampleRate;
paramFilt.BW = 1e3;
paramFilt.freq_central = f_Q;
paramFilt.order = 2;
paramFilt.gain = 1;
paramFilt.plot_flag = false;

FilteredToneQ = custom_filter(PD_Signal, paramFilt);
Spectrum_FilteredToneQ = abs(fftshift(fft(FilteredToneQ)));

figure,plot(time,FilteredToneQ)
figure,plot(frequency,10*log10(Spectrum_FilteredToneQ),'k')
title('Spectrum of the Filtered Pilot tone Q'), axis([0 1e4 -60 Inf])

%% SWEEP
% 
% 
% 
