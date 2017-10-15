%% Sweep P
SweepPts = 300;
SweepRangeMax = 3*Vpi.P;

% Minimum Transmitance
bias.I = 1.0 * Vpi.I; 
bias.Q = 1.0 * Vpi.Q;
bias.P = 0:SweepRangeMax/(SweepPts-1):SweepRangeMax; % Sweep
tone.P = 0;
E_OUT_I = mzm1(E_IN, 0, 0, bias.I, Vpi.I);
E_OUT_Q = mzm1(E_IN, 0, 0, bias.Q, Vpi.Q);
E_OUT_IQ = E_OUT_I + quad(E_OUT_Q, bias.P, tone.P);
PD.R = 0.95;     % Responsivity
PD_Power = PD.R*abs(E_OUT_IQ).^2;
figure('Position', [488 280 560 208]),
plot(PD_Power)
title('Sweep Bias P - I/Q Min'), ylabel('PD Power'), xlabel('Sample')
%%
% Quadrature Point
bias.I = 0.5 * Vpi.I; 
bias.Q = 0.5 * Vpi.Q;
bias.P = 0:SweepRangeMax/(SweepPts-1):SweepRangeMax; % Sweep
tone.P = 0;
E_OUT_I = mzm1(E_IN, 0, 0, bias.I, Vpi.I);
E_OUT_Q = mzm1(E_IN, 0, 0, bias.Q, Vpi.Q);
E_OUT_IQ = E_OUT_I + quad(E_OUT_Q, bias.P, tone.P);
PD.R = 0.95;     % Responsivity
PD_Power = PD.R*abs(E_OUT_IQ).^2;
figure('Position', [488 280 560 208]),
plot(PD_Power)
title('Sweep Bias P - I/Q Quad'), ylabel('PD Power'), xlabel('Sample')

% Maximum Transmitance
bias.I = 0.0 * Vpi.I; 
bias.Q = 0.0 * Vpi.Q;
bias.P = 0:SweepRangeMax/(SweepPts-1):SweepRangeMax; % Sweep
tone.P = 0;
E_OUT_I = mzm1(E_IN, 0, 0, bias.I, Vpi.I);
E_OUT_Q = mzm1(E_IN, 0, 0, bias.Q, Vpi.Q);
E_OUT_IQ = E_OUT_I + quad(E_OUT_Q, bias.P, tone.P);
PD.R = 0.95;     % Responsivity
PD_Power = PD.R*abs(E_OUT_IQ).^2;
figure('Position', [488 280 560 208]),
plot(PD_Power)
title('Sweep Bias P - I/Q Max'), ylabel('PD Power'), xlabel('Sample')