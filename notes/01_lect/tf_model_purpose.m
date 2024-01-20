% This is the worked out example given in class: transfer functions are
% similar, but their closed loop responses are not.

close all; clear variables;
format compact

% G = \frac{50}{(10s + 1)} \frac{10}{s^2 + 2s + 10}
G = tf(50, [10, 1]) * tf(10, [1, 2, 10]);

% G_1 = \frac{50}{(10s + 1)} \frac{1}{(0.3s+1)} \frac{1}{(0.1s+1}
G_1 = tf(50, [10, 1]) * tf(1, [0.3, 1]) * tf(1, [0.1, 1]);

% G_2 = \frac{10}{(2s+1)} \frac{10}{(s^2 + 2s + 10)}
G_2 = tf(10, [2, 1]) * tf(10, [1, 2, 10]);

% controller
C = 0.01 * tf([10, 1], [0.1, 0]);

L     = G*C;
L_1   = G_1*C;
L_2   = G_2*C;
L_cl  = feedback(L, 1);
L_1cl = feedback(L_1, 1);
L_2cl = feedback(L_2, 1);

figure();
bode(G, "b", G_1, "r--", G_2, "c--", logspace(-1,1.4,201)); grid("on")
legend("G", "G_1", "G_2", "Location", "northeast"); legend("boxoff")
title("open loop Bode diagram"); % saveas(gcf, "bode_ol.png")
figure();
step(L_cl, "b", L_1cl, "r--", L_2cl, "c--", 6); grid("on")
legend("G", "G_1", "G_2", "Location", "northeast"); legend("boxoff"); grid("on")
ylimit()
title("closed loop step response"); % saveas(gcf, "step_cl.png")
figure();
bode(L_cl, "b", L_1cl, "r--", L_2cl, "c--", logspace(-1,1.4,201)); grid("on")
legend("G", "G_1", "G_2", "Location", "southwest"); legend("boxoff")
title("closed loop Bode diagram"); % saveas(gcf, "bode_cl.png")

figure();
step(G, G_1, G_2, 10); legend("G", "G_1", "G_2", "Location", "northwest"); legend("boxoff"); grid("on")
title("open loop step response")
figure();
nyquist(L, L_1, L_2); legend("L", "L_1", "L_2", "Location", "northwest"); legend("boxoff")
title("(open loop) Nyquist plot")
