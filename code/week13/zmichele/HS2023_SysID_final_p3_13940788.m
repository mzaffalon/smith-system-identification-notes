function [p31_Tur_hat, p31_Tyr_hat, p31_G_hat, p31_omega, p31_r1, p31_r2, p32_theta_hat] = HS2023_SysID_final_p3_13940788()
%% Solution template for Problem 3

name = mfilename;
legi_num = str2double(name(end-7:end));

DEBUG = false;
%legi_num = 13940788;

%% Part 1:

%Define a reference input
u = idinput(255, "prbs");
p31_r1 = [u(end-3:end); repmat(u,4,1)];
p31_r2 = p31_r1;
assert(length(p31_r1) == 1024);

% Excite system and obtain data for part 1 and 2
[p31_u2, p31_y1, p32_r, p32_y]  = HS2023_SysID_final_p3_GenerateData(legi_num, p31_r1, p31_r2);

p31_r1 = p31_r1(5:end);
p31_y1 = p31_y1(5:end);
p31_r2 = p31_r2(5:end);
p31_u2 = p31_u2(5:end);
Ts = 1/255;

%%
% For part 1, do not use p32_r, p32_y

Tyr = fft(sum(reshape(p31_y1,255,4), 2)) ./ fft(sum(reshape(p31_r1,255,4), 2));
Tur = fft(sum(reshape(p31_u2,255,4), 2)) ./ fft(sum(reshape(p31_r2,255,4), 2));

p31_Tur_hat = Tur(2:128); % vector 127x1
p31_Tyr_hat = Tyr(2:128); % vector 127x1
p31_G_hat   = p31_Tyr_hat ./ p31_Tur_hat; % vector 127x1
p31_omega   = 2*pi*(1:127)'; % vector 127x1

if DEBUG
    figure(); t=tiledlayout(2,1);
    nexttile(); semilogy(p31_omega/(2*pi), abs(p31_G_hat));
    xlabel("freq [1/s]"); ylabel("|G|"); grid("on")
    nexttile(); plot(p31_omega/(2*pi), rad2deg(unwrap(angle(p31_G_hat))));
    xlabel("freq [rad/s]"); ylabel("angle [deg]"); grid("on")
    title(t, "Ghat")
end

disp("===== 3.1.a) =====");
disp("The scaled PRBS signals that I have chosen have white dense spectrum (persistently exciting) and small DC content. They are scaled to be in the range [-1,1] so as to maximize the signal power. Since we need to find the frequency response in the range [1,...,127]Hz, I use a periodic signal with period 255=127*2+1 to increase SNR.")

disp("===== 3.1.b) =====");
disp("I first average the responses Y1 and U2 over the 4 periods and then use ETFE to estimate Tyr=fft(Y1)./fft(R1) and Tur=fft(U2)./fft(R2). The joint input-output method relies on the two estimates Tyr and Tur to have different noise realizations.")

disp("===== 3.1.c) =====");
disp("Tyr and Tur are unbiased since they are ETFE, if it was not for the unknown transient; they are however asymptotically unbiased. G=Tyr/Tur is biased (not asymptotically unbiased) since the ratio does not preserve Gaussian statistics but it is consistent.")

disp("===== 3.1.e) =====");
disp("The unknown transient generates a biased estimate. I assume the two experiments cannot be run back to back. I would use the first experiment to determine the length of the transient using an impulse and the second to estimate G after discarding the transient. Since now I know the length of the transient, I would work in time-domain and estimate the frequency response at the frequencies 1..127Hz.");

%% Part 2:
% For part 2, do not use p31_u2, p31_y1 p31_r1, p31_r2

C = tf([0.5,0], [1,-0.9], Ts);
u = p32_r - lsim(C, p32_y);

Phiyu = [-toeplitz(p32_y(1:end-1), [p32_y(1),0,0]), ...
    toeplitz(u(1:end-1), [u(1),0])];

disp("===== 3.2.a) =====");
disp("ARX forms a linear regressor and gives unbiased estimates.");

disp("===== 3.2.b) =====");
disp("Since G is estimated with an ARX model, the estimate is unbiased, not only consistent, provided C(z) is correct.")

p32_theta_hat = Phiyu \ p32_y(2:end);

end
