% Malte Mueller 05.04.2019

%% Parameter

U = 1;
T = 5;
tau = 1;
dt = tau/100;

%% Aufgabe 1

% Berechnung der Lade- und Entladekurve

t1 = 0:dt:tau;              % Ladezeit
uy1 = U * (1-exp(-t1/tau)); % Ladefunktion

t2 = tau:dt:(tau+5*tau);    % Entladezeit
uy2 = U * (1-exp(-tau/tau)) * exp(-(t2-tau)/tau);   % Entladefunktion

uy = [uy1, uy2];    % Zusammensetzen der Lade- und Entladewerte
t = [t1, t2];

plot(t, uy); grid on; grid minor; %axis tight; 

% Beschriftung des Graphen

title('Ausgangsspannung u_y(t)');
xlabel('t');
ylabel('u_y(t)');

%% Aufgabe 2

dt = tau/100;
% Berechnung
for m = 1:10
    tau = m*tau/5;
    
    t1 = 0:dt:tau;
    uy1 = U * (1-exp(-t1/tau));
    
    t2 = tau:dt:((tau+5*tau));
    uy2 = U * (1-exp(-tau/tau)) * exp(-(t2-tau)/tau);
    
    uy = [uy1 uy2];
    t = [t1 t2];
    
    if (m == 1 || m == 3 || m == 10)
        plot(t, uy);
        hold on;
    end
    
    %legend(strcat('m = ', num2str(m)));
    pause(1.0);
end

% Beschriftung des Plots
title('Ausgangsspannung uy(t)');
legend('m = 1', 'm = 3', 'm = 10');
xlabel('t');
ylabel('uy');

%% Aufgabe 3

% Einstellen der Parameter
tau = tau;
t1 = linspace(0, tau, 1000);
t2 = linspace(tau, 4*T, 4000);
F = 1/T;
t = [t1, t2];
% figure;
% Approximation
for K = 1:100
    uyt = 0;

    for k = -K:K
        Hkf = 1 / (1 + 1i*2*pi*k*F*tau);
        Xk = exp(-1i*2*pi*k*F*(tau/2)) * (tau/T) * sinc(k*F*tau);
        % Xk = (U/T) * ((exp(-1i*2*pi*k*F*tau) - 1)/(-1i*2*pi*k*F));
        Yk = Hkf * Xk;
        uyt = uyt + Yk * exp(1i*2*pi*k*F*t);
    end
    if (K == 2 || K == 8 || K == 40 || K == 100)
        plot(t, uyt);
        hold on;
    end
    %legend(strcat('K =', num2str(K)));
    %pause(0.5);
end
% Beschriftung des Plots

title('Approximierte Ausgangsspannung uy(t)');
xlabel('t');
ylabel('u_y');
legend('K = 2', 'K = 8', 'K = 40', 'K = 100');

hold off;