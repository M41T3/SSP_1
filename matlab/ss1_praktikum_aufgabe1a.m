% Malte M?ller 05.04.2019

%% Parameters

U = 1;
T = 5;
T1 = 1;

%% Aufgabe 1 (c)FERK

% Einstellen der Parameter
tau = T1;
dt = T1/100;
% Berechnung
t1 = 0:dt:T1;
uy1 = U * (1-exp(-t1/tau));

t2 = T1:dt:(T1+5*tau);
uy2 = U * (1-exp(-T1/tau)) * exp(-(t2-T1)/tau);

uy = [uy1, uy2];
t = [t1, t2];

plot(t, uy);
% Beschriftung des Plots
title('Ausgangsspannung uy(t)');
xlabel('t');
ylabel('uy');

%% Aufgabe 2 (c)FERK

dt = T1/100;
% Berechnung
for m = 1:10
    tau = m*T1/5;
    
    t1 = 0:dt:T1;
    uy1 = U * (1-exp(-t1/tau));
    
    t2 = T1:dt:((T1+5*tau));
    uy2 = U * (1-exp(-T1/tau)) * exp(-(t2-T1)/tau);
    
    uy = [uy1 uy2];
    t = [t1 t2];
    
    if (m == 1 || m == 3 || m == 10)
        plot(t, uy);
        hold on;
    end
    
    %legend(strcat('m = ', num2str(m)));
    %pause;
end

% Beschriftung des Plots
title('Ausgangsspannung uy(t)');
legend('m = 1', 'm = 3', 'm = 10');
xlabel('t');
ylabel('uy');

%% Aufgabe 3 (c)FERK

% Einstellen der Parameter
tau = T1;
t1 = linspace(0, T1, 1000);
t2 = linspace(T1, T, 1000);
F = 1/T;
t = [t1 t2];
figure;
% Approximation
for K = 1:100
    uyt = 0;
    for k = -K:K
        Hkf = 1 / (1 + 1i*2*pi*k*F*tau);
        Xk = exp(-1i*2*pi*k*F*(T1/2)) * (T1/T) * sinc(pi*k*F*T1);
        Yk = Hkf * Xk;
        uyt = uyt + Yk * exp(1i*2*pi*k*F*t);
    end
    if (K == 2 || K == 8 || K == 40 || K == 100)
        plot(t, real(uyt));
        hold on;
    end
    %legend(strcat('K =', num2str(K)));
    %pause;
end
% Beschriftung des Plots

title('Approximierte Ausgangsspannung uy(t)');
xlabel('t');
ylabel('uy');
legend('K = 2', 'K = 8', 'K = 40', 'K = 100');

hold off;