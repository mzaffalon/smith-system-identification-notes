close all; clear variables;
format compact

%%

A = tf([1,-0.5,0.36], [1,0,0], 1);
B = tf([0.8,-0.65], [1,0,0], 1);
C = tf([1,1], [1,0], 1);
Niter = 2000;
sigma = 0.5;
G = minreal(B/A);
H = minreal(C/A);
u = prbs(5,2*127)';
N = length(u);
y_true = lsim(G, u);
invC = minreal(1/C);
Cm1 = minreal(C-1);

%%

thetas_e = zeros(Niter, 4);
thetas_Le = zeros(Niter, 4);

for i = 1:Niter
    e = sigma*randn(N,1);
    v = lsim(H, e);
    y = y_true + v;
    yF = lsim(invC, y);
    uF = lsim(invC, u);

    Phiyu = [-toeplitz(yF(1:end-1), [yF(1); 0]), ...
        toeplitz(uF(1:end-1), [uF(1); 0])];

    Cm1yF = lsim(Cm1, yF);
    theta = Phiyu \ (y(2:end) - Cm1yF(2:end));
    thetas_e(i,:) = theta';
    theta = Phiyu \ yF(2:end);
    thetas_Le(i,:) = theta';
end
disp(mean(thetas_Le))
disp(std(thetas_Le))
disp(mean(thetas_e))
disp(std(thetas_e))
