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
Vpi.I = 3.75;              % Half-wave voltage, Typ: 3.75 Vpp
Vpi.Q = 3.75;
Vpi.P = 3.75;
IL = 0;                 % Insertion Loss,  Typ: 8.5 dB
ER_I = 10^( 30 /10);    % Extinction Ratio, Typ: 30 dB
ER_Q = 10^( 30 /10);

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

figure,
subplot(211),plot(time,Signal_I)
title('Signal I'), ylabel('Samples'), axis([0 1e-4 -1.1 1.1])
subplot(212),plot(time,Signal_Q)
title('Signal Q'), ylabel('Samples'), axis([0 1e-4 -1.1 1.1])

figure('Position',[618 412 300 280]),
plot(Signal_I, Signal_Q, 'o-' ), xlabel('Signal I'), ylabel('Signal Q')
title('QPSK Generated Signal (Normalized)')

%% BW limitation
% Transfer function for the filter response:
%
% $$ TF_{BW}(f) = \exp(-\log(\sqrt{2})(\frac{f}{0.75*2})^2) $$
%
% (Not used)
if false
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

Vpp_I = 0.2; % tone amplitude (peak-to-peak)
Vpp_Q = 0.2;
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

%% Photodetector - Average optical power
PD.Var = 10^(-60/10); % -50 dB of Noise floor
PD.R = 0.95;     % Responsivity
PD_Power = PD.R*abs(E_OUT_IQ).^2;
PD_Noise_floor = PD.Var*randn(size(E_OUT_IQ));

%% Heater Transfer Function
paramFilt.SampleRate = SampleRate;
paramFilt.BW = 1e4;
paramFilt.freq_central = 0;
paramFilt.order = 1;
paramFilt.gain = 1;
paramFilt.plot_flag = false;
PD_Signal = custom_filter(PD_Power, paramFilt) + PD_Noise_floor;
Spectrum_PD_Signal = abs(fftshift(fft(PD_Signal)));
figure('Position',[384 442 800 320]),
subplot(121),plot(time,PD_Signal)
title('Photodectector Received Signal'), xlabel('Samples')
subplot(122),plot(frequency,10*log10(Spectrum_PD_Signal),'k')
title('Photodectector Received Signal'), axis([0 2e4 -60 Inf])
xlabel('Frequency'), ylabel('Power Spectral Density')

%% Pilot Tone Filters
% Ideal case, where all the Bias points are optimized.
plot_pilot_tones

%% SWEEP
%
SweepPts = 300;
SweepRangeMax = 3*Vpi.I;
