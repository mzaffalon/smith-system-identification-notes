% The ARX regressor is formed from the noisy measured data y(k), not from
% the true noiseless data. Is the estimate biased? Here I plot

% Model: y(k) = B(z)/A(z)*u(k) + 1/A(z)*e(k)
% One-step ahead predictor: yhat = (1-A(z))y(k) + B(z)u(k)

close all; clear variables;
format compact

%% this is the example given in class

%rng(14);
b = 0.5;
a = -0.7;
G = tf([b], [1, a], 1);
sigma = 0.25;
Niter = 1000;
u = 2*prbs(7,2*127)'-1;
N = length(u);
y_true = lsim(G,u);
Phi_true = make_Phi(y_true,1,u,1);

%% H1 = 1: statistics

H1 = tf(1, 1, 1);
thetas_fin = zeros(Niter,4);
for i = 1:Niter
    e = sigma * randn(N,1);
    y = y_true + lsim(H1, e);
    Phi = make_Phi(y,1,u,1);
    theta_tn = Phi_true \ y(2:end);
    theta_nn = Phi \ y(2:end);
    thetas_fin(i,:) = [theta_tn', theta_nn'];
end

disp("means")
disp(mean(thetas_fin))
disp("std")
disp(std(thetas_fin))

%% with H2

H2 = tf(1, [1, a], 1);
thetas_fin = zeros(Niter,4);
for i = 1:Niter
    e = sigma * randn(N,1);
    y = y_true + lsim(H2, e);
    Phi = make_Phi(y,1,u,1);
    theta_tn = Phi_true \ y(2:end);
    theta_nn = Phi \ y(2:end);
    thetas_fin(i,:) = [theta_tn', theta_nn'];
end

disp("means")
disp(mean(thetas_fin))
disp("std")
disp(std(thetas_fin))

%%

rng(14);
e = sigma * randn(N,1);
y = y_true + lsim(H2, e);
Phi = make_Phi(y,1,u,1);

thetas_inc = zeros(N-1,4);
for i = 10:N-1
    theta_tn = Phi_true(1:i,:) \ y(2:i+1);
    theta_nn = Phi(1:i,:) \ y(2:i+1);
    thetas_inc(i,:) = [theta_tn', theta_nn'];
end

figure(); hold("on"); plot(0:N-1, y_true); plot(0:N-1, y); grid("on")
legend("$y_{true}$", "$y_{noisy}$", "Interpreter","latex"); legend("boxoff"); xlim([0,N])

figure(); tiledlayout(2,1)
nexttile; hold("on"); plot(10:N-1, thetas_inc(10:end,2:2:end), "-", "LineWidth",1.5); plot(10:N-1, b*ones(N-10,1))
legend("$\hat{\theta}_{tn}$", "$\hat{\theta}_{nn}$", "b", ...
    "Location","southeast", "Interpreter","latex")
legend("boxoff"); title("b"); xlim([0,N])

nexttile; hold("on"); plot(10:N-1, thetas_inc(10:end,1:2:end), "-", "LineWidth",1.5); plot(10:N-1, a*ones(N-10,1)); 
legend("$\hat{\theta}_{tn}$", "$\hat{\theta}_{nn}$", "a", ...
    "Location","southeast", "Interpreter","latex")
legend("boxoff"); title("a"); xlim([0,N]); xlabel("length of regressor")

%% private functions

function Phi = make_Phi(y,n,u,m)
Phi = [toeplitz(-y(1:end-1), [-y(1); zeros(n-1,1)]), ...
    toeplitz(u(1:end-1), [u(1); zeros(m-1,1)])];
end
