%% Plot Quad IQ, Out of Phase
% Quadrature 
bias.I = 0.5 * Vpi.I;
bias.Q = 0.5 * Vpi.Q;
bias.P = 0.9 * Vpi.P; 

E_OUT_I = mzm1(E_IN, RF.I, tone.I, bias.I, Vpi.I);
E_OUT_Q = mzm1(E_IN, RF.Q, tone.Q, bias.Q, Vpi.Q);
E_OUT_IQ = E_OUT_I + quad(E_OUT_Q, bias.P, tone.P);

figure('Position',[618 412 300 280]),
plot(E_OUT_IQ, '.'), xlabel('Real part'), ylabel('Imaginary part')
title('QPSK Optical Signal - I min'), axis([-1.5 1.5 -1.5 1.5])

PD_Power = PD.R*abs(E_OUT_IQ).^2;
PD_Noise_floor = PD.Var*randn(size(E_OUT_IQ));
paramFilt.SampleRate = SampleRate;
paramFilt.BW = 1e4;
paramFilt.freq_central = 0;
paramFilt.order = 1;
paramFilt.gain = 1;
paramFilt.plot_flag = false;
PD_Signal = custom_filter(PD_Power, paramFilt) + PD_Noise_floor;

plot_pilot_tones


%%


