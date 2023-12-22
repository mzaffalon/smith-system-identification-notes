function [p21_thetaHat1, p21_thetaHat2, p21_thetaHat3, p21_bias1, ...
    p21_covMat1, p21_MSE1, p21_thetaHat1_reg, p21_thetaHat2_reg, ...
    p21_thetaHat3_reg, p22_ghat_ls_100, p22_n_cv, p22_ghat_ls_n_cv, p22_ghat_tc_cv] ...
    = HS2023_SysID_final_p2_13940788()
%
%% Outputs:
% Part 1 outputs:
%   p21_thetaHat1:    plant parameters of system 1 (10 x 1)
%   p21_thetaHat2:    plant parameters of system 2 (10 x 1)
%   p21_thetaHat3:    plant parameters of system 3 (10 x 1)
%   p21_bias1:        bias vector of system 1 (10 x 1)
%   p21_covMat1:      covariance matrix of system 1 (10 x 10)
%   p21_MSE1:         mean squared error of system 1 (1 x 1)
%   p21_thetaHat1_reg: regularized joint estimates for systems 1 (10 x 1)
%   p21_thetaHat2_reg: regularized joint estimates for systems 1 (10 x 1)
%   p21_thetaHat3_reg: regularized joint estimates for systems 1 (10 x 1)
%
% Part 2 outputs:
%  p22_ghat_ls_100:   best weighted least-squares estimate at order n=100
%                     (given all three datasets) (100 x 1)
%  p22_n_cv:          best length of weighted least-squares estimates
%                     from orders n=1:100 (through hold-out cross-validation)
%  p22_ghat_ls_n_cv:  best least-squares estimate from orders n=1:100
%                     (through hold-out cross-validation) (p22_n_cv x 1)
%  p22_ghat_tc_cv:       best TC-regularized least-squares estimate at order
%                     n=100 (through hold-out cross-validation) (100 x 1)
%

%% Begin routines
% Extract Legi from Filename
name=mfilename;
LegiNumber= name(end-7:end);

DEBUG = false;
%LegiNumber = 13940788;

[Y1, Y2, Y3, dataSet1, dataSet2, dataSet3, var_ratio] = HS2023_SysID_final_p2_GenerateData(LegiNumber);
sigma_ratio = sqrt(var_ratio);

%% Part 1

var_e = 0.5^2;
[p21_thetaHat1,p21_thetaHat1_reg,p21_bias1,p21_covMat1,p21_MSE1] = dumb_LS_response(Y1, var_e);
[p21_thetaHat2,p21_thetaHat2_reg,~,~,~] = dumb_LS_response(Y2, var_e);
[p21_thetaHat3,p21_thetaHat3_reg,~,~,~] = dumb_LS_response(Y3, var_e);

disp("The James-Stein estimator is guaranteed to provide a biased minimum MSE estimate for the problem Y=theta_0+e when N>2 by an appropriate shrinkage of the solution towards the origin.")

%% Part 2

order = 100;
[Phi1, y1] = makePhiy(dataSet1, order, order);
[Phi2, y2] = makePhiy(dataSet2, order, order);
[Phi3, y3] = makePhiy(dataSet3, order, order);
p22_ghat_ls_100 = [Phi1; Phi2/sigma_ratio; Phi3] \ [y1; y2/sigma_ratio; y3];

if DEBUG
    figure(); hold("on"); stem(Phi1 \ y1); stem(Phi2 \ y2);
    stem(Phi3 \ y3);
    stem(p22_ghat_ls_100, "LineWidth",2);
end

%%

taumax = 100;
M = cell(taumax,1);
E = zeros(taumax,1);
for order = 1:taumax
    [Phi1, y1]       = makePhiy(dataSet1, order, order);
    [Phi2, y2]       = makePhiy(dataSet2, order, order);
    [Phi_val, Y_val] = makePhiy(dataSet3, order, order);

    %[Phi1, y1]       = makePhiy(dataSet1, taumax, order);
    %[Phi2, y2]       = makePhiy(dataSet2, taumax, order);
    %[Phi_val, Y_val] = makePhiy(dataSet3, taumax, order);

    theta = [Phi1; Phi2/sigma_ratio] \ [y1; y2/sigma_ratio];
    M{order} = theta;
    E(order) = mean((Y_val - Phi_val*theta).^2);
end

figure(221); plot(1:taumax, E); xlabel("order"); ylabel("mean squared error")
title("Mean squared error on the validation data")

[~,p22_n_cv] = min(E);
p22_ghat_ls_n_cv = M{p22_n_cv};

%%

n = 100;
[Phi1, y1]       = makePhiy(dataSet1, n, n);
[Phi2, y2]       = makePhiy(dataSet2, n, n);
[Phi_val, Y_val] = makePhiy(dataSet3, n, n);
Phi_id           = [Phi1; Phi2/sigma_ratio];
Y_id             = [y1;   y2/sigma_ratio];

alphas = logspace(-1,0,41);
alphas = alphas(1:end-1);
gammas = logspace(0,2,31);
%MM = cell(length(gammas), length(alphas));
EE = zeros(length(gammas), length(alphas));
GA = cell(length(gammas), length(alphas));
for ia = 1:length(alphas)
    alpha = alphas(ia);
    invP = TCKernel_inv(alpha, n);
    for ig = 1:length(gammas)
        gamma = gammas(ig);
        theta = (transpose(Phi_id)*Phi_id + gamma*invP) \ (transpose(Phi_id)*Y_id);
        %MM{ig, ia} = theta;
        EE(ig, ia) = mean((Y_val - Phi_val*theta).^2);
        GA{ig, ia} = [gamma, alpha];
    end
end

[~,II] = min(EE, [], "all"); % I cannot figure out
%p22_ghat_tc_cv = MM{II};
gamma_opt = GA{II}(1);
alpha_opt = GA{II}(2);
invP = TCKernel_inv(alpha_opt, n);
Phi = [Phi_id; Phi_val];
Y   = [Y_id; Y_val];
p22_ghat_tc_cv = (transpose(Phi)*Phi + gamma_opt*invP) \ (transpose(Phi)*Y);

if DEBUG
    figure(); contourf(alphas, gammas, E); colorbar(); xlabel("alpha"); ylabel("gamma")
    figure(); surf(alphas, gammas, E); colorbar(); xlabel("alpha"); ylabel("gamma")
    figure(); hold("on"); stem(p22_ghat_ls_n_cv), stem(p22_ghat_tc_cv); legend("LS", "reg")
end

%% All done!

end

%% Student-defined functions

function [theta,theta_reg,bias,cov,MSE] = dumb_LS_response(Y, var_e)
theta = Y;
r = 8*var_e / norm(Y,2)^2; % hard-coded number of elements N-2=8
theta_reg = (1-r)*theta;
bias = zeros(10,1);
cov = var_e*eye(10);
MSE = 10*var_e;
end


function [Phi, y] = makePhiy(data, taumax, order)
u = data(:,1);
% For y discard the first taumax elements. NaNs in the first taumax
% elements count towards the transient, so that we remove only the NaNs
% after the transient.
y = data(taumax+1:end,2);
% find and remove the NaNs
m = find(~isnan(y));
y = y(m);
% construct the toeplitz matrix; the initial conditions are not known, so
% start from the taumax+1 entry to compute g(i), i=0..order
Phi = toeplitz(u(taumax+1:end), flip(u(taumax-order+1:taumax+1)));
Phi = Phi(m,2:end);
end


function invP = TCKernel_inv(alpha, n)
diag0 = alpha .^ (-(1:n)') / (1-alpha);
diag0(2:n-1) = diag0(2:n-1) * (1+alpha);
diag1 = -alpha .^ (-(1:n-1)') / (1-alpha);
invP = diag(diag0, 0) + diag(diag1, +1) + diag(diag1, -1);
end
