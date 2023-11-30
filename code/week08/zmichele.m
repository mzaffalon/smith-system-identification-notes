close all; clear variables;
format compact

%%

N = 100; % period length in points
tau_max = 80;
d = 40; % number of sinusoidals
G = tf([11/10, 5/6], [1, -5/8, 0], 1);

u = sinusoidal_input(N, 2, d);
g0 = single_shot(G, u, tau_max, 0); % noiseless

sigma = sqrt(0.05);

%%

M = 200;
gest_avg = zeros(tau_max+1,M);
% accumulate the averages
gest_avg(:,1) = single_shot(G, u, tau_max, sigma);
for i = 2:M
    gest_avg(:,i) = gest_avg(:,i-1) + single_shot(G, u, tau_max, sigma);
end
gest_avg = gest_avg ./ (1:M);
e_avg = transpose(vecnorm(gest_avg - g0, 2));

figure(); 
loglog(1:M, e_avg); hold("on"); loglog(1:M, e_avg(1) ./ sqrt(1:M)); loglog(1:M, e_avg(1) ./ (1:M))
grid("on"); xlabel("number of averages"); ylabel("error"); legend("error", "1/sqrtM", "1/M"); legend("boxoff")
title("Averaging method (error)")

%%

m = 30;
gest_N = zeros(tau_max+1,m-1);
gest_Nr = zeros(tau_max+1,m-1);
for i = 1:m-1
    w = sinusoidal_input(N, i+1, d);
    gest_N(:,i) = single_shot(G, w, tau_max, sigma);
    wr = random_input((i+1)*tau_max);
    gest_Nr(:,i) = single_shot(G, wr, tau_max, sigma);
end
e_N  = transpose(vecnorm(gest_N  - g0, 2));
e_Nr = transpose(vecnorm(gest_Nr - g0, 2));

figure();
loglog(2:m, e_N); hold("on"); loglog(2:m, e_Nr);
loglog(2:m, sqrt(2)*e_N(1) ./ sqrt(2:m))
grid("on"); xlabel("number of periods"); ylabel("error"); 
legend("periodic input", "rand input", "1/sqrtM"); legend("boxoff")
title("m periods (error)")


%% user defined functions

function [u, freq] = sinusoidal_input(N, m, d)
% generate the input signal consisting of d sinusoidals
freq = 2*pi/N*round(N/2*(1:d)'/(d+1)); % equally spaced frequencies excluding extremes
thetas = pi*(1:d)'.*((1:d)'-1)/d; % Schroeder phase
t = (0:N-1)';
u = zeros(N,1);
for p = 1:d
    u = u + sin(freq(p)*t + thetas(p));
end
% normalize the input signal to fall in the range [-2,2]
u = 2*u/norm(u,Inf);
u = repmat(u,m,1);
end


% This seems to fix only when all frequencies are included
function [u, freq] = sinusoidal_input2(N, m, d)
[u, freq] = idinput([N 1 m], "sine", [0 1], [-2 2], [d 10 1]);
end


function u = random_input(N)
    u = rand(N,1) - 0.5; % remove constant term
    u = 2*u/norm(u,"inf");
end

%%%%%%%%%%
function gest = single_shot(G, u, tau_max, sigma)
y = lsim(G, u) + sigma * randn(length(u),1);
Phiu = toeplitz(u, [u(1) zeros(1, tau_max)]);
gest = Phiu \ y;
%figure(); plot(y-Phiu*gest); grid("on"); xlabel("points"); ylabel("diff"); title("y-Phi_u*gest")
end
