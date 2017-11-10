%% Bessel Functions - Test
Pin = 1;
Vpi=3.75;

%% Sweep V_DC
Vdc=0:(3*Vpi)/99:(3*Vpi);

Pout = (Pin/2) .* (1+cos(pi*Vdc/Vpi));

figure,
plot(Vdc, Pout)
xlabel('V_{dc}'), ylabel('P_{out}/P_{in}')
title('Sweep V_{DC}')

xticks([0 Vpi 2*Vpi 3*Vpi]);
xticklabels({'0','V_{\pi}','2V_{\pi}','3V_{\pi}'});
yticks([0 0.25 0.5 0.75 1]);

%% Plot Fundamental mode of Pilot Tone
A = .2;      % ~200mV
f=2e3; w=2*pi*f;

t=0:4*pi*(1/w)/99:4*pi*(1/w);
Vac = A*sin(w*t);

figure,
plot(t,Vac)
title('Fundamental mode of Pilot Tone')
xlabel('t'), ylabel('Amplitude')
xlabel('t')

%% Plot 2nd Harmonic of Pilot Tone
f=4e3; w=2*pi*f; 

t=0:8*pi*(1/w)/99:8*pi*(1/w);
Vac = A*sin(w*t);

figure,
plot(t,Vac)
title('2nd Harmonic of Pilot Tone')
xlabel('t')

%% Quadrature Bias - Fundamental
Vac = .0533*Vpi;    % ~200mV
f = 2e3; w = 2*pi*f;
t = 0:2*pi*(1/w)/99:2*pi*(1/w);
Vdc = (Vpi/2)*ones(1,100);

Vin = Vdc + Vac*sin(w*t);
Pout = Pin*cos(pi*Vin/Vpi);

figure,
plot(Vin, Pout)
xlabel('V_{in}'), ylabel('P_{out}/P_{in}')
title('Quadrature Bias point - Fundamental')

xticks([Vpi/2-0.1 Vpi/2 Vpi/2+0.1]);
xticklabels({'V_{\pi}/2 - 0.1','V_{\pi}/2','V_{\pi}/2 + 0.1'});
yticks([-0.1 0 0.1]);

%% Minimum Bias point - 2nd Harmonic
Vdc = (Vpi)*ones(1,100);
Vin = Vdc + Vac*sin(w*t);
Pout = Pin*cos(pi*Vin/Vpi);

figure,
plot(Vin, Pout)
xlabel('V_{in}'), ylabel('P_{out}/P_{in}')
title('Minimum Bias point - 2nd Harmonic')

xticks([Vpi-.0999 Vpi Vpi+.1]);
xticklabels({'V_{\pi}- 0.1','V_{\pi}','V_{\pi}+ 0.1'});
yticks([-1 -.998 -.996]);

%% Maximum Bias point - 2nd Harmonic
Vdc = zeros(1,100);
Vin = Vdc + Vac*sin(w*t);
Pout = Pin*cos(pi*Vin/Vpi);

figure,
plot(Vin, Pout)
xlabel('V_{in}'), ylabel('P_{out}/P_{in}')
title('Maximum Bias point - 2nd Harmonic')

xticks([-.1 0 .1]);
xticklabels({'-0.1','0','0.1'});
yticks([.996 .998 1]);

%% Harmonic Analysis
Pin = 1;    % 1mW
Vpi = 3.75; % 3.75 V
f = 2e3;    % 2 kHz
w = 2*pi*f;

k = 20 *(2*pi);      % number of periods
% it determines the speed of the Vdc sweep
% higher number means slower sweep

T = k*(1/f);
fs = 10000; % Sample rate = 100 kHz
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

Pout0 = (Pin/2)*cos(pi*Vdc/Vpi).*J0;
Pout1 = Pin*sin(pi*Vdc/Vpi).*J1.*sin(w*t);
Pout2 = Pin*cos(pi*Vdc/Vpi).*J2.*cos(2*w*t);
Pout3 = Pin*sin(pi*Vdc/Vpi).*J3.*sin(3*w*t);
Pout4 = Pin*cos(pi*Vdc/Vpi).*J4.*cos(4*w*t);

Poutb = Pin/2 + Pout0 - Pout1 + Pout2 - Pout3 + Pout4;

figure,
plot(Vdc, Pout)
ylabel('P_{out}'), xlabel('V_{DC}')
title(['Harmonic Analysis, V_{AC}= ' num2str(Vac/Vpi,3) ' V_{\pi}= ' num2str(Vac,3) 'V'])
xticks([0 Vpi/2 Vpi 3*Vpi/2 2*Vpi]);
xticklabels({'0','V_{\pi}/2','V_{\pi}','3V_{\pi}/2','2V_{\pi}'});

figure,
hold on
plot(t,J0)
plot(t,J1)
plot(t,J2)
plot(t,J3)
plot(t,J4)
legend('J_0(x)','J_1(x)','J_2(x)','J_3(x)','J_4(x)')
xlabel('t'), title(['Bessel Functions, x= V_{AC}= ' num2str(Vac,3) 'V'])
xticks([0 T/2 T]);
xticklabels({'0','T/2','T'});

figure,
hold on
plot(t,Pout0)
plot(t,-Pout1)
plot(t,Pout2)
plot(t,-Pout3)
plot(t,Pout4)
legend('P_{out,0}','-P_{out,1}','P_{out,2}','-P_{out,3}','P_{out,4}')
xlabel('t'), title('Harmonic Analysis - Bessel Functions')

%%
figure,
subplot(211)
hold on
plot(Vdc, Pout)
plot(Vdc, Poutb)
ylabel('Bessel P_{out}'), xlabel('V_{DC}')
title('Bessel Functions Approximation')
legend('P_{out}','Bessel 4th order approx')
xticks([0 Vpi/2 Vpi 3*Vpi/2 2*Vpi]);
xticklabels({'0','V_{\pi}/2','V_{\pi}','3V_{\pi}/2','2V_{\pi}'});
subplot(212)
error=zeros(1,length(Vdc));
for i=1:length(Vdc)
    if Pout(i) > 0
        error(i) = abs(Pout(i)-Poutb(i))./Pout(i);
    else
        error(i)=nan;
    end
end
plot(Vdc, error)
xticks([0 Vpi/2 Vpi 3*Vpi/2 2*Vpi]);
xticklabels({'0','V_{\pi}/2','V_{\pi}','3V_{\pi}/2','2V_{\pi}'});

%%
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
xlim([0 2*Vpi]), ylim([-120 0]), grid
xlabel('V_{DC}'), ylabel('P_{avg} (dBm)')
title('Detected DC and Harmonic Contents')
xticks([0 Vpi/2 Vpi 3*Vpi/2 2*Vpi]);
xticklabels({'0','V_{\pi}/2','V_{\pi}','3V_{\pi}/2','2V_{\pi}'});

%% Harmonics Analysis
% Harmonic distorsion caused by pilot tone's amplitude applied 
% Vdc = Vpi/2 (quadrature points)
Vac = .0533*Vpi;     % ~200mV
Pout = (Pin/2)*( 1 + cos( pi/2 + (pi*Vac/Vpi).*sin(w*t) ) );
L=length(t);
ff = (100)*(0:L/2-1)/(2*pi);
S = fft(Pout);
mag = abs(S(1:L/2)/L);
figure,
plot(ff, mag)
axis([0 1e4 0 .2])
ylabel('|S_{out}|'), xlabel('f')
title('Power density spectrum, V_{AC}= 200mV')

Vac = 0.5*Vpi;
Pout = (Pin/2)*( 1 + cos( pi/2 + (pi*Vac/Vpi).*sin(w*t) ) );
% L=length(t);
% ff = (fs/2)*(0:L/2-1)/(2*pi);
S = fft(Pout);
mag = abs(S(1:L/2)/L);
figure, plot(ff, mag)
axis([0 1e4 0 .25])
ylabel('|S_{out}|'), xlabel('f')
title(['Power density spectrum, V_{AC}= 0.5V_{\pi}= ' num2str(.5*Vpi,3) 'V'])

Vac = 1*Vpi;
Pout = (Pin/2)*( 1 + cos( pi/2 + (pi*Vac/Vpi).*sin(w*t) ) );
% L=length(t);
% ff = (fs/2)*(0:L/2-1)/(2*pi);
S = fft(Pout);
mag = abs(S(1:L/2)/L);
figure, plot(ff, mag)
axis([0 1e4 0 .25])
ylabel('|S_{out}|'), xlabel('f')
title('Power density spectrum, V_{AC}= V_{\pi}')

%% Harmonic Analysis
% Harmonic distorsion caused by pilot tone's amplitude applied 
% Vdc = 0 V (maximum transmission point)
Vac = 2*.0533*Vpi;     % ~200mV
Pout = (Pin/2)*( 1 + cos( (pi*Vac/Vpi).*sin(w*t) ) );
% L=length(t);
% ff = (fs/2)*(0:L/2-1)/(2*pi);
S = fft(Pout);
mag = abs(S(1:L/2)/L);
figure,
plot(ff, mag)
axis([0 1e4 0 .025])
ylabel('|S_{out}|'), xlabel('f')
title(['Power density spectrum, V_{AC}= ' num2str(Vac,3) 'V'])

Vac = .5*Vpi;
Pout = (Pin/2)*( 1 + cos( (pi*Vac/Vpi).*sin(w*t) ) );
S = fft(Pout);
mag = abs(S(1:L/2)/L);
figure, plot(ff, mag)
axis([0 1e4 0 .25])
ylabel('|S_{out}|'), xlabel('f')
title(['Power density spectrum, V_{AC}= 0.5V_{\pi}= ' num2str(Vac,3) 'V'])

Vac = Vpi;
Pout = (Pin/2)*( 1 + cos( (pi*Vac/Vpi).*sin(w*t) ) );
S = fft(Pout);
mag = abs(S(1:L/2)/L);
figure, plot(ff, mag)
axis([0 1e4 0 .25])
ylabel('|S_{out}|'), xlabel('f')
title('Power density spectrum, V_{AC}= V_{\pi}')
%%
%% Color map
Pin = 1;    % 1mW
Vpi = 3.75; % 3.75 V
f = 2e3;    % 2 kHz
w = 2*pi*f;

k = 2 *(2*pi);      % number of periods
% it determines the speed of the Vdc sweep
% higher number means slower sweep

T = k*(1/f);
fs = 1000; % Sample rate = 1 kHz (for Color map)
t=0:T/(fs-1):T;

% Vac = .0267*Vpi;     % ~100mV
% Vac = .0533*Vpi;     % ~200mV
% Vac = .1067*Vpi;     % ~400mV
% Vac = .1867*Vpi;     % ~700mV
% Vac = .2667*Vpi;     % ~1V
% Vac = .3733*Vpi;     % ~1.4V
Vac = Vpi; % 3.75V

Vdc = 0:2*Vpi/(fs-1):2*Vpi;

figure,
[X,Y] = meshgrid(t,Vdc);
Pout = (Pin/2)*( 1 + cos( pi*Y/Vpi + (pi*Vac/Vpi).*sin(w*X) ) );
surf(X,Y, Pout,'LineStyle','none')
xlabel('t'), ylabel('V_{DC}'), zlabel('P_{out}')
title(['Harmonic Analysis, V_{AC}= ' num2str(Vac/Vpi,3) 'V_{\pi}= ' num2str(Vac,3) 'V'])

yticks([0 Vpi/2 Vpi 3*Vpi/2 2*Vpi]);
yticklabels({'0','V_{\pi}/2','V_{\pi}','3V_{\pi}/2','2V_{\pi}'});

%%
%% Plot Bessel Functions
w=1; t=0:4*pi*(1/w)/999:4*pi*(1/w);
Vac = 0.1;
phi = Vac*ones(1,length(t));
J0 = besselj(0,phi);
J1 = besselj(1,phi);
J2 = besselj(2,phi);
J3 = besselj(3,phi);
J4 = besselj(4,phi);
figure,
hold on
plot(t,J0)
plot(t,J1)
plot(t,J2)
plot(t,J3)
plot(t,J4)
legend('J_0(x)','J_1(x)','J_2(x)','J_3(x)','J_4(x)')
title('Bessel Functions, x= 0.1')
xticks(); xticklabels({});

%% Plot Bessel Functions
w=1; t=0:4*pi*(1/w)/999:4*pi*(1/w);
Vac = 0.1*Vpi;
phi = (Vac)*ones(1,length(t));
J0 = besselj(0,phi);
J1 = besselj(1,phi);
J2 = besselj(2,phi);
J3 = besselj(3,phi);
J4 = besselj(4,phi);
figure,
hold on
plot(t,J0)
plot(t,J1)
plot(t,J2)
plot(t,J3)
plot(t,J4)
legend('J_0(x)','J_1(x)','J_2(x)','J_3(x)','J_4(x)')
title(['Bessel Functions, x= ' num2str(Vac/Vpi,3) ' V_{\pi}'])
xticks(); xticklabels({});

%% Plot Bessel Functions
w=1; t=0:1/999:1;
Vac = 1;
phi = (Vac)*ones(1,length(t));
J0 = besselj(0,phi);
J1 = besselj(1,phi);
J2 = besselj(2,phi);
J3 = besselj(3,phi);
J4 = besselj(4,phi);
figure,
hold on
plot(t,J0)
plot(t,J1)
plot(t,J2)
plot(t,J3)
plot(t,J4)
legend('J_0(x)','J_1(x)','J_2(x)','J_3(x)','J_4(x)')
title('Bessel Functions, x= 1')
xticks(); xticklabels({});

%% Plot Bessel Functions
x=0:4*pi/999:4*pi;

J0 = besselj(0,x);
J1 = besselj(1,x);
J2 = besselj(2,x);
J3 = besselj(3,x);
J4 = besselj(4,x);

figure,
hold on
plot(x,J0)
plot(x,J1)
plot(x,J2)
plot(x,J3)
plot(x,J4)
legend('J_0(x)','J_1(x)','J_2(x)','J_3(x)','J_4(x)')
xlabel('x'), title('Bessel Functions')

xticks([0 pi 2*pi 3*pi 4*pi]);
xticklabels({'0','\pi','2\pi','3\pi','4\pi'});
