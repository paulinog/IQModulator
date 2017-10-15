%% Plot Sweep Phase
% Optimum bias point for phase is 0.5Vpi, which means in-quadrature point
%
% Pilot Tones were not applied
bias.I = 1.0 * Vpi.I; % Optimized for QPSK Modulation
bias.Q = 1.0 * Vpi.Q;
E_OUT_I = mzm1(E_IN, RF.I, 0, bias.I, Vpi.I);
E_OUT_Q = mzm1(E_IN, RF.Q, 0, bias.Q, Vpi.Q);
bias.P = 0.5 * Vpi.P;
E_OUT_IQ = E_OUT_I + quad(E_OUT_Q, bias.P, 0);
bias.P = 0.3 * Vpi.P;
E_OUT_IQ1 = E_OUT_I + quad(E_OUT_Q, bias.P, 0);
bias.P = 0.4 * Vpi.P;
E_OUT_IQ2 = E_OUT_I + quad(E_OUT_Q, bias.P, 0);
bias.P = 0.6 * Vpi.P;
E_OUT_IQ3 = E_OUT_I + quad(E_OUT_Q, bias.P, 0);
bias.P = 0.7 * Vpi.P;
E_OUT_IQ4 = E_OUT_I + quad(E_OUT_Q, bias.P, 0);

figure('Position',[618 412 480 280]),
hold on
plot(E_OUT_IQ, '.')
plot(E_OUT_IQ2, '.')
plot(E_OUT_IQ1, '.')
hold off
xlabel('Real part'), ylabel('Imaginary part')
title('QPSK Optical Signal - Model 1'), axis([-2 2 -2 2])
legend('Optimum, 0.5','0.4','0.3','Location','bestoutside')

figure('Position',[618 412 480 280]),
hold on
plot(E_OUT_IQ, '.')
plot(E_OUT_IQ3, '.')
plot(E_OUT_IQ4, '.')
hold off
xlabel('Real part'), ylabel('Imaginary part')
title('QPSK Optical Signal - Model 1'), axis([-2 2 -2 2])
legend('Optimum, 0.5','0.6','0.7','Location','bestoutside')