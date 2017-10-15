%% Sweep I
% In this case, P (phase of IQ quadrature) does not have any influence in 
% the optical output power variance.
SweepPts = 300;
SweepRangeMax = 3*Vpi.I;
% Bias Voltage Point
bias.I = 0:SweepRangeMax/(SweepPts-1):SweepRangeMax; % Sweep
bias.Q = 1.0 * Vpi.Q * ones(1,length(bias.I));
bias.P = 0.5 * Vpi.P; % Optimum point
tone.P = 0;
% QPSK Modulation, Model 1
E_OUT_I = mzm1(E_IN, 0, 0, bias.I, Vpi.I);
E_OUT_Q = mzm1(E_IN, 0, 0, bias.Q, Vpi.Q);
E_OUT_IQ = E_OUT_I + quad(E_OUT_Q, bias.P, tone.P);
PD.R = 0.95;     % Responsivity
PD_Power = PD.R*abs(E_OUT_IQ).^2;
figure('Position', [488 280 560 208]),
plot(PD_Power)
title('Sweep Bias I'), ylabel('PD Power'), xlabel('Sample')
