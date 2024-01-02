close all; clear variables;
format compact

%%
M      = 63;
m      = 128;
m_trans = 3;
N     = M * m;
sigma = 0.2;
rng(15);

Ts    = 1.0;
G     = tf([0.1, 0], [1, -2.2, 2.42, -1.87, 0.7225], Ts);
H     = tf([0.5, -0.5*0.9], [1, -0.25], Ts);

w     = idinput(M, "rgs");
w     = idinput(M, "prbs");
u     = [zeros(m_trans*M,1); repmat(w, m, 1)];
u_per = repmat(w,m+m_trans,1);

%%
e     = sigma * randn(N+m_trans*M,1);
v     = lsim(H, e);
t     = (-m_trans*M:N-1)';
% t is irrelevant for the simulation because the system is an LTI
y     = lsim(G, u)     + v;
y_per = lsim(G, u_per) + v;

%%

idx = m_trans*M + (-1*M:2*M)';
tt = t(idx);
figure(); tlay = tiledlayout(3,1, "TileSpacing","tight");
nexttile(); plot(tt, u(idx), "o-");            title("u(k)");      xlim([min(tt),max(tt)])
nexttile(); plot(tt, y(idx), "o-");            title("y(k)");      xlim([min(tt),max(tt)])
nexttile(); plot(tt, y_per(idx)-y(idx), "o-"); title("r(k)");      xlim([min(tt),max(tt)])
title(tlay, "input and output signals")

saveas(gcf, "transient-noise-time-domain.png")

%%
f = (0:M)'/M/Ts;
idx = (f > 0 & f <= 0.5); % skip w=0
f = f(idx);
omega = 2*pi*f;
Gtrue = squeeze(freqresp(G, omega));
ftw = fft(w);
ftw = ftw(idx);
m_plots = 2.^(0:log2(m));
E = zeros(m,1);

figure(); tlay = tiledlayout(2,1);
nexttile(1); semilogy(f, abs(Gtrue), "LineWidth",2, "DisplayName","true G");
hold("on"); grid("on"); %legend("boxoff"); legend("Location","west")
nexttile(2); plot(f, rad2deg(unwrap(angle(Gtrue))), "LineWidth",2);
hold("on"); grid("on"); ylim([-600, 100])

for i = 1:m
    rz = reshape(y(m_trans*M+1:m_trans*M+i*M),[],i);
    ftz = fft(mean(rz,2));
    Ghat = ftz(idx) ./ ftw;
    E(i) = norm(Ghat-Gtrue,2);
    if any(i == m_plots)
        nexttile(1); semilogy(f, abs(Ghat), "DisplayName",sprintf("Ghat m=%d",i))
        nexttile(2); plot(f, rad2deg(unwrap(angle(Ghat))))
    end
end
nexttile(1); xlabel("freq [1/s]"); ylabel("|Gest|")
nexttile(2); xlabel("freq [1/s]"); ylabel("angle [deg]")

figure(); plot(f, abs(ftw)); grid("on");
xlabel("freq [1/s]"); ylabel("|U(e^{j\omega_n})|"); title("PSD of input signal")

figure(); loglog(1:m, E, "LineWidth",2); hold("on");
loglog(1:m, E(1)./(1:m)'); loglog(1:m, sqrt(m)*E(end)./sqrt(1:m)');
xlabel("m"); ylabel("error"); legend("error", "1/m", "1/sqrt(m)"); legend("boxoff")
title("||Ghat - Gtrue||_2"); grid("on"); ylim([0.9*min(E), 1.1*max(E)])
saveas(gcf, "transient-noise-convergence.png")

%% Average the experiment Niter times

Niter = 100;
y_true = lsim(G, u);
E = zeros(m,1);
for k = 1:Niter
    e = sigma * randn(N+m_trans*M,1);
    v = lsim(H, e);
    y_n = y_true + v;
    for i = 1:m
        rz = mean(reshape(y_n(m_trans*M+1:m_trans*M+i*M),[],i),2);
        ftz = fft(rz);
        Ghat = ftz(idx) ./ ftw;
        E(i) = E(i) + norm(Ghat-Gtrue,2);
    end
end

E = E / Niter;
figure(); loglog(1:m, E, "LineWidth",2); hold("on");
loglog(1:m, E(1)./(1:m)'); loglog(1:m, sqrt(m)*E(end)./sqrt(1:m)');
xlabel("m"); ylabel("error"); legend("error", "1/m", "1/sqrt(m)"); legend("boxoff")
title("||Ghat - Gtrue||_2"); grid("on"); ylim([0.9*min(E), 1.1*max(E)])

saveas(gcf, "transient-noise-convergence-avg.png")
