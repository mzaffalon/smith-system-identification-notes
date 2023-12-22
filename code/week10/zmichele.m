close all; clear variables;
format compact

%%

LegiNumber = 13940788;
[u, y] = HS2023_SysID_Exercise_08_GenerateData(LegiNumber);

theta_LS = solve_regLS(u, y, 20, 0, 0.5);
theta_TC = solve_regLS(u, y, 20, 100, 0.6);

%%

N_id = length(u) * 0.7;
u_id  = u(1:N_id);
y_id  = y(1:N_id);
u_val = u(N_id+1:end);
y_val = y(N_id+1:end);

rs = (2:30)';
alphas = (0:0.05:0.95)';

% It is not necessary to save all the thetas: we only need to find the
% minimum and it could be done in one pass.
M = cell(length(rs), length(alphas));
E = zeros(length(rs), length(alphas));
for ir = 1:length(rs)
    for ia = 1:length(alphas)
        M{ir, ia} = solve_regLS(u_id, y_id, rs(ir), 50, alphas(ia));
        E(ir, ia) = residuals_regLS(M{ir, ia}, u_val, y_val, rs(ir), 50, alphas(ia));
    end
end

figure(); contourf(alphas, rs, E); colorbar(); xlabel("alpha"); ylabel("order")
figure(); surf(alphas, rs, E); colorbar(); xlabel("alpha"); ylabel("order")

%%

% find the minimum
[~,I] = min(E, [], "all");
stem_fig = figure(); stem(theta_LS); hold("on"); stem(theta_TC); stem(M{I})
xlabel("parameter number"); ylabel("value")
legend("LS", "regul. LS (TC)", "cross-validation"); legend("boxoff")

disp("Best alpha = " + alphas(floor(I / length(rs)) + 1))
disp("Best r     = " + rs(mod(I, length(rs))))

%% test correctness of inverse
n = 12;
alpha = 0.2;
assert(norm(TCKernel_inv(alpha, n) * TCKernel(alpha, n) - eye(n), Inf) < 1e-12)
assert(norm(TCKernel(alpha, n) * TCKernel_inv(alpha, n) - eye(n), Inf) < 1e-12)

%% local function
function invP = TCKernel_inv(alpha, n)
diag0 = alpha .^ (-(1:n)') / (1-alpha);
diag0(2:n-1) = diag0(2:n-1) * (1+alpha);
diag1 = -alpha .^ (-(1:n-1)') / (1-alpha);
invP = diag(diag0, 0) + diag(diag1, +1) + diag(diag1, -1);
end

function theta = solve_regLS(u, y, r, gamma, alpha)
Phi = toeplitz(u, [u(1) zeros(1, r-1)]);
if gamma == 0
    theta = Phi \ y;
    return
end
if alpha == 0
    % invP has Inf on the diagonals
    theta = zeros(r,1);
    return
end
invP = TCKernel_inv(alpha, r);
% generalised SVD is about 100 times slower than the other two solutions
% [U,~,X,C,~] = gsvd(Phi, sqrt(gamma) * chol(invP));
% theta = transpose(X) \ (transpose(U*C)*y);

% solve the problem [Phi, chol(invP)] theta - [y, 0]: this is the fastest
% solution!
theta = [Phi; sqrt(gamma) * chol(invP)] \ [y; zeros(r,1)];
% avoid computing the inverse, but not really better...
% theta = (transpose(Phi)*Phi + gamma*invP) \ (transpose(Phi)*y);
end

function res = residuals_regLS(theta, u, y, r, gamma, alpha)
Phi = toeplitz(u, [u(1) zeros(1, r-1)]);
res = norm(y - Phi*theta, 2)^2;
if alpha ~= 0 && gamma ~= 0
    invP = TCKernel_inv(alpha, r);
    res = res + gamma * transpose(theta) * invP * theta;
end
% the official solution does not normalize
%res = sqrt(res / length(y));
end


%% used only for testing the inverse
function P = TCKernel(alpha, n)
P = zeros(n,n);
alpha_n = 1;
for i = 1:n
    alpha_n = alpha_n * alpha;
    P(i,1:i) = alpha_n;
    P(1:i-1,i) = alpha_n;
end
end
