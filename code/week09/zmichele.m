close all; clear variables;
format compact

%%
N = 200;
sigma_u = 1;
sigma_e = 0.1;

% The model is y = B/A*u + 1/A*e

G = tf([0.6, 0.3, -0.05, 0], [1, 0.25, -0.2, 0.1, 0.05], 1);
H = tf(1, G.Denominator{1}, G.Ts);

u_id = sigma_u * randn(N,1);
y_id_true = lsim(G, u_id);
v_id = lsim(H, sigma_e * randn(N,1));
y_id_noise = y_id_true + v_id;

%%

p_fr = figure();
bode(G); grid("on");
title("transfer function and estimates")
l_ft = legend("G");

n_coeff = 12;
err_id  = zeros(n_coeff,1);
err_val = zeros(n_coeff,1);
u_val = sigma_u * randn(N,1);
y_val_noise = lsim(G, u_val) + lsim(H, sigma_e * randn(N,1));

for i = 1:n_coeff
    % estimate the polynomials
    n_a = i; n_b = i;
    Phi = [-toeplitz(y_id_noise(1:end-1), [y_id_noise(1) zeros(1,n_a-1)]), ...
        toeplitz(u_id(1:end-1), [u_id(1) zeros(1,n_b-1)])];
    theta = Phi \ y_id_noise(2:end);

    %disp(norm(y(2:end) - Phi*theta,2))
    b = transpose(theta(n_a+1:end));
    a = transpose([1; theta(1:n_a)]);
    Gest = tf(b, a, G.Ts);
    y_est = lsim(Gest, u_id);
    % figure(); plot(y_est); grid("on")
    err_id(i) = norm(y_est-y_id_noise, 2)^2/N;

    y_val_est = lsim(Gest, u_val);
    err_val(i) = norm(y_val_est-y_val_noise, 2)^2/N;

    if n_a > 4
        hold("on")
        bode(Gest)
        %legend(l_fr, "est ")
    end
end
figure()
semilogy(1:n_coeff, err_id); hold("on")
semilogy(1:n_coeff, err_val); grid("on")
xlabel("number of coefficients"); ylabel("error")
legend("ident", "valid"); legend("boxoff")
title("convergence error")
