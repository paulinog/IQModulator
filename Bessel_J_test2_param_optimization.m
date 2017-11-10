%% Parameters Optimization
Pin = 1;    % 1mW
Vpi = 3.75; % 3.75 V
f = 2e3;    % 1 kHz
w = 2*pi*f;

k = 20 *(2*pi);      % number of periods
% it determines the speed of the Vdc sweep
% higher number means slower sweep

T = k*(1/f);
fs = 100000; % Sample rate = 1000*f/k
t=0:T/(fs-1):T;

% Vac = .0267*Vpi;     % ~100mV
Vac = .0533*Vpi;     % ~200mV
% Vac = .1067*Vpi;     % ~400mV
% Vac = .1867*Vpi;     % ~700mV
% Vac = .2667*Vpi;     % ~1V
% Vac = .3733*Vpi;     % ~1.4V

Vdc = 0:2*Vpi/(fs-1):2*Vpi;

Pout = (Pin/2)*( 1 + cos( pi*Vdc/Vpi + (pi*Vac/Vpi).*sin(w*t) ) );

phi = (Vac)*ones(1,length(t));
J0 = besselj(0,phi);
J1 = besselj(1,phi);
J2 = besselj(2,phi);
J3 = besselj(3,phi);
J4 = besselj(4,phi);

figure,
plot(Vdc, Pout)
ylabel('P_{out}'), xlabel('V_{DC}')
title(['Harmonic Analysis, V_{AC}= ' num2str(Vac/Vpi,3) ' V_{\pi}= ' num2str(Vac,3) 'V'])
xticks([0 Vpi/2 Vpi 3*Vpi/2 2*Vpi])
xticklabels({'0','V_{\pi}/2','V_{\pi}','3V_{\pi}/2','2V_{\pi}'})

figure,
hold on
plot(t,J0)
plot(t,J1)
plot(t,J2)
plot(t,J3)
plot(t,J4)
legend('J_0(x)','J_1(x)','J_2(x)','J_3(x)','J_4(x)')
xlabel('t'), title(['Bessel Functions, x= V_{AC}= ' num2str(Vac,3) 'V'])
xticks([0 T/2 T])
xticklabels({'0','T/2','T'})

Pavg0 = ((Pin/2) + (Pin/2)*cos(pi*Vdc/Vpi).*J0)*T;
Pavg1 = Pin*sin(pi*Vdc/Vpi).*J1.*(1-cos(w*T))/(w*T);
Pavg2 = Pin*cos(pi*Vdc/Vpi).*J2.*sin(2*w*T)/(2*w*T);
Pavg3 = Pin*sin(pi*Vdc/Vpi).*J3.*(1-cos(3*w*T))/(3*w*T);
Pavg4 = Pin*cos(pi*Vdc/Vpi).*J4.*sin(4*w*T)/(4*w*T);

figure,
hold on
plot(Vdc, 10*log10(abs(Pavg0)))
plot(Vdc, 10*log10(abs(Pavg1)))
plot(Vdc, 10*log10(abs(Pavg2)))
plot(Vdc, 10*log10(abs(Pavg3)))
plot(Vdc, 10*log10(abs(Pavg4)))
legend('P_{out,0}','P_{out,1}','P_{out,2}','P_{out,3}','P_{out,4}','Location','eastoutside')
xlabel('V_{DC}'), ylabel('P_{avg} (dBm)')
title('Detected DC and Harmonic Contents')
xticks([0 Vpi/2 Vpi 3*Vpi/2 2*Vpi])
xticklabels({'0','V_{\pi}/2','V_{\pi}','3V_{\pi}/2','2V_{\pi}'})
