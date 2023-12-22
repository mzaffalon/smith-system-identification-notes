function [p1_theta1_a, p1_theta1_b,p1_theta2, p1_y2_hat] = HS2023_SysID_final_p1_13940788()

%% Solution template for Final Problem 1

name = mfilename;
LegiNumber = str2double(name(end-7:end));

DEBUG = false;
%LegiNumber = 13940788;

[p1_u1, p1_y1, p1_u2, p1_y2, p1_w] = HS2023_SysID_final_p1_GenerateData(LegiNumber);

%% Write your solutions here

%% Task 1
Phiyu = [-toeplitz(p1_y1(1:end-1), [p1_y1(1); 0]), ...
    toeplitz(p1_u1(1:end-1), [p1_u1(1); 0])];
p1_theta1_a = Phiyu \ p1_y1(2:end);

C = tf([1,1], [1,0], 1);
invC = minreal(1/C);
yC = lsim(invC, p1_y1);
uC = lsim(invC, p1_u1);
PhiyuC = [-toeplitz(yC(1:end-1), [yC(1); 0]), ...
    toeplitz(uC(1:end-1), [uC(1); 0])];
p1_theta1_b = PhiyuC \ yC(2:end);

if DEBUG
    yC_hat = PhiyuC * p1_theta1_b;
    yhat = lsim(C, yC_hat);
    figure(); hold("on"); plot(p1_y1(2:end)); plot(Phiyu * p1_theta1_a); plot(yhat);
    grid("on"); legend("true", "p1\_theta\_a", "p1\_theta\_b"); legend("boxoff")
    figure(); hold("on"); plot(p1_y1(2:end) - Phiyu * p1_theta1_a); plot(p1_y1(2:end) - yhat);
    grid("on"); legend("true-p1\_theta\_a", "true-p1\_theta\_b"); legend("boxoff")
end

disp("===== 1 =====");
disp("(c) The first problem is a standard ARX and the least-squares estimate is known to be unbiased and minimal variance (amongst all that are unbiased) as seen in class. The second problem can be reduced to an ARX after letting yF = C^{-1}y and uF=C^{-1}u where C(z)=1-z^{-1}. Strictly speaking, the transformation is not justified because C(z) has a zero on the unit circle and is therefore not stably invertible. A numerical test however shows that it still works.")
disp("(d) An ARMAX model structure solved with least-squares gives an estimate that has minimum variance but is biased. This is because ARMAX regressor is pseudolinear: the unknown theta enters also the regressor.")


%% Task2

% implement yH_hat = (1-D)*yH + F*uH + H2*wH with yH = H_1^{-1}y
H1 = tf([1,1,0], [1, p1_theta1_a(1:2)'], 1);
invH1 = minreal(1/H1);
y2H = lsim(invH1, p1_y2);
u2H = lsim(invH1, p1_u2);
wH  = lsim(invH1, p1_w);

PhiyuwH = [-toeplitz(y2H(1:end-1), [y2H(1); 0]), ...
    toeplitz(u2H(1:end-1), [u2H(1); 0]), ...
    wH(1:end-1)];
p1_theta2 = PhiyuwH \ (y2H(2:end) - wH(2:end));  % [d1,d2,f1,f2,h2]

D  = tf([1, p1_theta2(1:2)'], [1,0,0], 1);
F  = tf(p1_theta2(3:4)', [1,0,0], 1);
H2 = tf([1,p1_theta2(5),0], [1,0,0], 1);
p1_y2_hat = lsim(F, u2H) + lsim(H2, wH) + lsim(minreal(H1 - D), y2H);

if DEBUG
    H1yH = lsim(minreal(H1-1), y2H);
    theta = PhiyuwH \ (p1_y2(2:end) - wH(2:end) - H1yH(2:end));
    if norm(theta - p1_theta2, Inf) > 1e-10
        disp("solutions not compatible")
    end
end

disp("===== 1 =====");
disp("(c) With the transformation yH = H^{-1}y, uH=H^{-1}u, wH=H^{-1}w, we have an ARX. wH behaves like a second control input with a pass-through term that accounts for an offset in the one-step ahead predictor: yhat(k|k-1) = [(1-D)y_H + Fu_H + z^{-1}wH(k)] + wH(k) and G_2=F/D.")
disp("(d) With the transformation above, the system is an ARX with unbiased estimate (so more than asymptotically unbiased) and minimum variance.")

figure(111)
plot(p1_y2); hold("on"); plot(p1_y2_hat)
legend("signal p1\_y2", "one-step ahead p1\_y2\_hat"); legend("boxoff"); grid("on")
xlabel("time"); ylabel("signal"); title("Comparison signal vs one-step ahead estimate")

end
