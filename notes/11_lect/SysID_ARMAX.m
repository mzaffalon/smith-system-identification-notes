% The name ARMAX was already taken...

close all; clear variables;
format compact

%%

Ts = 1;
Bnum = [1, 0.5];       B = tf(Bnum, [1, 0, 0], Ts); % no constant term (no feed-through)
Anum = [1, -1.7, 0.9]; A = tf(Anum, [1, 0, 0], Ts); % monic
Cnum = [1, -3, 1];     C = tf(Cnum, [1, 0, 0], Ts); % monic
G = minreal(B/A);
H = minreal(C/A);

sigma = 1;
Nexps = 128;
O = 5;
K = 31;
O = 7;
K = 127;
lambda = 0.0;

u = 2*prbs(O,K) - 1;
ytrue = lsim(G, u);

%% reproducibility

% What is the influence of the solver? Starting with the same (u,y),
% compute a solution from different starting conditions.

Xs = zeros(Nexps, 6); % collect the coefficients

e = sigma * randn(K,1);
v = lsim(H, e);
y = ytrue + v;

for i = 1:Nexps
    x0 = [zeros(6,1); randn(K,1)]; % change only the starting condition
    x = ARMAX(u,y,x0,lambda);
    Xs(i,:) = x(1:6);
end

figure(); boxplot(Xs, "Labels", ["b1", "b2", "a1", "a2", "c1", "c2"])
hold("on"); plot(1:6, [Bnum, Anum(2:end), Cnum(2:end)], "o")
title(sprintf("ARMAX coefficients: solver reproducibility (PRBS O=%d, K=%d)", O, K))
figure(); boxplot(Xs(:,1:4), "Labels", ["b1", "b2", "a1", "a2"])
hold("on"); plot(1:4, [Bnum, Anum(2:end)], "o")
title(sprintf("ARMAX coefficients: solver reproducibility (PRBS O=%d, K=%d)", O, K))

%%

Xs = zeros(Nexps, 6); % collect the coefficients
Es = zeros(Nexps, 2); % collect the residulas: we expect \epsilon = e

for i = 1:Nexps
    e = sigma * randn(K, 1);
    v = lsim(H, e);
    y = ytrue + v;
    [x, fval] = ARMAX(u,y,zeros(6+K,1),lambda);
    Xs(i,:) = x(1:6);
    Es(i,:) = [dot(e(7:end),e(7:end))+lambda*dot(e(1:6),e(1:6)), fval]/K; % true vs estimated residuals
end

figure(); boxplot(Xs, "Labels", ["b1", "b2", "a1", "a2", "c1", "c2"])
hold("on"); plot(1:6, [Bnum, Anum(2:end), Cnum(2:end)], "o")
title(sprintf("ARMAX coefficients (PRBS O=%d, K=%d)", O, K))

figure(); t = tiledlayout(2,1);
nexttile(); plot(Es); grid("on"); xlabel("experiment iteration"); ylabel("error")
legend("true", "estimated"); legend("Box","off");
err = Xs(:,5:6)-Cnum(2:3);
nexttile(); plot(err); grid("on"); xlabel("experiment iteration"); ylabel("estimated - true value")
title(t, sprintf("residuals (PRBS O=%d, K=%d)", O, K))

%%

Bhat = tf(    x(1:2)',  [1, 0],    Ts);
Ahat = tf([1, x(3:4)'], [1, 0, 0], Ts);
Chat = tf([1, x(5:6)'], [1, 0, 0], Ts);
Ghat = minreal(Bhat/Ahat);
Hhat = minreal(Chat/Ahat);
ehat = x(7:end);
vhat = lsim(Ghat, ehat);
yhat = lsim(Ghat, u) + vhat; % should I add the estimated noise too when I plot?

figure(); stem(ytrue); hold("on"); stem(yhat)
legend("y noisy", "yhat"); legend("boxoff")
title("Signal: noisy vs estimate")

figure(); stem(v); hold("on"); stem(vhat)
legend("v", "vhat"); legend("boxoff")
title("Noise: true vs estimate")

%% definition of the minimizer

function [x,fval] = ARMAX(u,y,x0,lambda)
PhiTyu = [toeplitz(u(1:end-1), [u(1); 0]), ...
    toeplitz(-y(1:end-1), [-y(1); 0])];

[x,fval] = fmincon(@(x) ARMAXobjective(x,lambda),x0, ...
    [],[],[],[],[],[],@(x)ARMAXconstraint(x,y,PhiTyu), ...
    optimoptions("fmincon", "Display","none", ...
    "SpecifyObjectiveGradient",false, ...
    "SpecifyConstraintGradient",false));
end


function [f,g] = ARMAXobjective(x,lambda) % x = [theta; e]
f = dot(x(7:end), x(7:end)) + lambda * dot(x(1:6),x(1:6));
if nargout == 2
    % attemp to use the gradients too, but the solver becomes _SLOWER_
    g = 2*[lambda*x(1:6); x(7:end)];
end
end


function [c,ceq,gc,gceq] = ARMAXconstraint(x,y,PhiTyu)
theta = x(1:6);
e = x(7:end);
PhiTe = toeplitz(e(1:end-1), [e(1); 0]);
ceq = y(2:end) - [PhiTyu, PhiTe] * theta - e(2:end);
c = [];
if nargout == 4
    Nd = length(e)-1;
    T = toeplitz(zeros(Nd,1), [0; 1; zeros(Nd-1,1)]);
    K = toeplitz([theta(5:6); zeros(Nd-2,1)], [theta(5); zeros(Nd,1)]);
    gceq = -[PhiTyu, PhiTe, K+T]';
    gc = [];
end
end
