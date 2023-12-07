% 

close all; clear variables;
format compact

%%
rng(10) % fix the RNG seed
order = 50;
G = tf([1,2,1], [1,-1.5,0.56,0],1);
g = impulse(G, order-1); % determine the true 50 elements of the impulse response
N = 1000;
sigma = 10;
u = randn(N,1);
e = sigma * randn(N,1);
y_true = lsim(G,u);
y = y_true + e;

figure(); plot(y_true); hold("on"); plot(y); grid("on")
legend("y true", "y noisy"); legend("boxoff")
xlabel("time"); ylabel("signal"); title("y true vs noisy")

%% split the data into identification and validation sets

N_id  = round(0.7 * N);
N_val = N - N_id;
u_id  = u(1:N_id);
y_id  = y(1:N_id);
u_val = u(N_id+1:end);
y_val = y(N_id+1:end);

thetaLS = solve_regLS(u, y, 50, 0);
figure(); stem(0:order-1, g); hold("on"); stem(0:order-1, thetaLS)
xlabel("coefficient index"); ylabel("amplitude"); title("true impulse response vs LS estimate")

%%

gammas = logspace(-1,2,201)';
Ngammas = length(gammas);
M_ridge = cell(Ngammas,1);
E = zeros(Ngammas,3);
BCE = zeros(Ngammas,3);
for i = 1:Ngammas
    invKernel = gammas(i)*eye(order);
    M_ridge{i} = solve_regLS(u_id, y_id, order, invKernel);
    % to choose the optimal gamma, we rely on the prediction error (the
    % residuals without the regularization term). We plot also the
    % regularization term and the full cost function for comparison.
    [res, prederr, regul] = residuals_regLS(M_ridge{i}, u_val, y_val, invKernel);
    E(i,:) = [res, prederr, regul] / (sigma^2*N_val);
    [b,c,m] = BCM(g, u_val, invKernel);
    BCE(i,:) = [b,c,m];
end

figure(); t = tiledlayout(2,1);
nexttile(); semilogx(gammas, E(:,2), "LineWidth",2); hold("on"); semilogx(gammas, E(:,[1,3]))
legend("prediction err", "cost func", "regul");
%nexttile(); semilogx(gammas, E(:,2))
%legend("prediction err");
legend("boxoff"); grid("on")
xlabel("gamma"); ylabel("residuals |Y-\Phi\theta|_2^2")

nexttile(); loglog(gammas, BCE); grid("on")
xlabel("gamma"); ylabel("contributions"); title("Regularized method: bias, cov and MSE")
legend("bias", "cov", "MSE", "Location", "south"); legend("boxoff")
title(t, sprintf("Ridge regression with %d points", N))

%%

% rs = (2:30)';
% 
% % It is not necessary to save all the thetas: we only need to find the
% % minimum and it could be done in one pass.
% M = cell(length(rs), length(gammas));
% E = zeros(length(rs), length(gammas));
% for ir = 1:length(rs)
%     for ia = 1:length(gammas)
%         M{ir, ia} = solve_regLS(u_id, y_id, rs(ir), 50, gammas(ia));
%         E(ir, ia) = residuals_regLS(M{ir, ia}, u_val, y_val, rs(ir), 50, gammas(ia));
%     end
% end
% 
% figure(); contourf(gammas, rs, E); colorbar(); xlabel("alpha"); ylabel("order")
% figure(); surf(gammas, rs, E); colorbar(); xlabel("alpha"); ylabel("order")
% 
% %%
% 
% % find the minimum
% [~,I] = min(E, [], "all");
% stem_fig = figure(); stem(theta_LS); hold("on"); stem(theta_TC); stem(M{I})
% xlabel("parameter number"); ylabel("value")
% legend("LS", "regul. LS (TC)", "cross-validation"); legend("boxoff")
% 
% disp("Best alpha = " + gammas(floor(I / length(rs)) + 1))
% disp("Best r     = " + rs(mod(I, length(rs))))


%% local function
function invP = TCKernel_inv(alpha, n)
diag0 = alpha .^ (-(1:n)') / (1-alpha);
diag0(2:n-1) = diag0(2:n-1) * (1+alpha);
diag1 = -alpha .^ (-(1:n-1)') / (1-alpha);
invP = diag(diag0, 0) + diag(diag1, +1) + diag(diag1, -1);
end

function theta = solve_regLS(u, y, r, invKernel)
Phi = toeplitz(u, [u(1); zeros(r-1,1)]);
if invKernel == 0
    theta = Phi \ y;
    return
end
theta = (transpose(Phi)*Phi + invKernel) \ (transpose(Phi)*y);
end

function [res, prederr, regul] = residuals_regLS(theta, u, y, invKernel)
r = length(theta);
Phi = toeplitz(u, [u(1) zeros(1, r-1)]);
prederr = norm(y - Phi*theta, 2)^2;
if invKernel == 0
    res = prederr;
    regul = 0;
    return
end
regul = transpose(theta) * invKernel * theta;
res = prederr + regul;
end

function [norm_bias, tr_cov, MSE] = BCM(theta0, u, invKernel)
r = length(theta0);
Phi = toeplitz(u, [u(1) zeros(1, r-1)]);
PhiTPhi = transpose(Phi)*Phi;
K = (PhiTPhi + invKernel) \ eye(r); % instead of using inv
norm_bias = norm(K*invKernel*theta0,2)^2;
tr_cov = trace(K*PhiTPhi*K);
MSE = norm_bias + tr_cov;
end
