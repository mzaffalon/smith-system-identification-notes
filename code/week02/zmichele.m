close all; clear variables;
format compact

%%
load("SysID_Exercise_1.mat");

mu = 0.5;
sigma = sqrt(0.3);

X = [x, 2*x.^2-1] / sigma;
b_ML = (y - mu) / sigma;
theta_ML = X\b_ML;
disp("theta_ML ="); disp(theta_ML)

mu_theta = [1.3, 0.9];
sigma_theta = sqrt(0.02);

% stick the a priori observation at the end of the measurement vector
X_MAP = [X; eye(2) / sigma_theta];
b_MAP = [b_ML; transpose(mu_theta) / sigma_theta];
theta_MAP = X_MAP \ b_MAP;
disp("theta_MAP ="); disp(theta_MAP)

% test on the validation data
X_v = [x_v, 2*x_v.^2-1];
E_ML = mean((X_v*theta_ML + mu - y_v).^2);
E_MAP = mean((X_v*theta_MAP + mu - y_v).^2);
disp("E_ML ="); disp(E_ML)
disp("E_MAP ="); disp(E_MAP)
