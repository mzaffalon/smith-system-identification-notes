close all; clear variables;
format compact

%%
M = 256;
m = 4;
N = m * M;
sigma = 0.0;
Ts = 1.0;

G = tf([0.1, 0], [1, -2.2, 2.42, -1.87, 0.7225], Ts);
H = tf([0.5, -0.5*0.9], [1, -0.25], Ts);

f = (0:M-1)'/M/Ts;
msk = (0 < f) & (f <= 0.5);
f = f(msk);
omega = 2*pi*f;
Gfresp = squeeze(freqresp(G, omega));
Hfresp = squeeze(freqresp(H, omega));

w = 2 * rand(M,1) - 1;
u = repmat(w,m,1);
e = sigma * randn(N,1);
y = lsim(G, u) + lsim(H, e);

% We want to take the mean _before_ the ratio
y_avg = mean(reshape(y, M, m), 2);
fty = fft(y_avg, [], 1);
ftw = fft(w); % there is input power only at the frequencies detemined by w
Gest = fty ./ ftw;
Gest = Gest(msk);
fprintf(1, "std(Gest-Gfresp)      = %f\n", std(Gest-Gfresp))

% using MATLAB etfe, gives the same result
Gest_etfe = etfe(iddata(y_avg, w, Ts));
Gest_etfe = squeeze(freqresp(Gest_etfe, omega));
fprintf(1, "std(Gest_etfe-Gfresp) = %f\n", std(Gest_etfe-Gfresp))

%%
figure()
tiledlayout(2,1);
nexttile
loglog(f, abs(Gfresp), "DisplayName", "|Gfresp|")
hold("on"); grid("on")
%loglog(f, abs(Hfresp), "DisplayName", "Hfresp")
loglog(f, abs(Gest),             "DisplayName", "|Gest|")
loglog(f, abs(Gest_etfe),        "DisplayName", "|Gest (efte)|")
loglog(f, abs(Gfresp-Gest),      "DisplayName", "|Gfresp-Gest|")
loglog(f, abs(Gfresp-Gest_etfe), "DisplayName", "|Gfresp-Gest\_efte|")
legend("Location", "west")
legend("boxoff")
xlabel("f [sec^{-1}]"); ylabel("magnitude")
xlim([0, 0.5]); ylim([5e-3, 3])

nexttile
semilogx(f, rad2deg(unwrap(angle(Gfresp))),    "DisplayName", "Gfresp")
hold("on"); grid("on")
semilogx(f, rad2deg(unwrap(angle(Gest))),      "DisplayName", "Gesp")
semilogy(f, rad2deg(unwrap(angle(Gest_etfe))), "DisplayName", "Gest (efte)")
legend("boxoff")
xlabel("f [sec^{-1}]"); ylabel("phase [deg]")
xlim([0, 0.5]); ylim([-600, 30])
