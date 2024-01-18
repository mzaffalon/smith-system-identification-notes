close all; clear variables;
format compact

%%

N  = 64;
sigma = 0.0;
rng(15);

Ts = 1.0;
% G  = tf([0.1, 0], [1, -2.2, 2.42, -1.87, 0.7225], Ts);
A = [2.2    -1.21    0.935  -0.7225;
     2        0        0        0;
     0        1        0        0;
     0        0      0.5        0];
B = [0.25; 0; 0; 0];
C = [0, 0, 0.2, 0];
D = 0;
G = tf(ss(A,B,C,D,Ts));
H = tf([0.5, -0.5*0.9], [1, -0.25], Ts);

% Compute the frequency response on a uniform grid on the entire unit
% circle, since we will need to extend it anyway.
w  = 2*pi*(0:N-1)'/(N*Ts);
G0 = squeeze(freqresp(G, w));

% add noise
e     = sigma * randn(2*N,1);
v     = lsim(H, e);
ftv   = fft(v(N+1:end));

Gest  = G0 + ftv;

%%

hhat = real(ifft(Gest));
q    = floor(N/2)+1;
Hhat = hankel(hhat(1:q), hhat(q:end));

[U,S,V] = svd(Hhat);
figure(); plot(diag(S), "o-"); grid("on")
n_xhat = 4; % we know the true system has n_x=4 from ss(G).A
S1 = S(1:n_xhat,1:n_xhat);
U1 = U(:,1:n_xhat);
V1 = V(:,1:n_xhat);
disp(norm(U1*S1*V1' - Hhat))

J3 = [1 zeros(1,q-1)];
Chat = J3 * U1; % = U1(1,:);
J1 = [eye(q-1), zeros(q-1,1)];
J2 = [zeros(q-1,1), eye(q-1)];
Ahat = (J1*U1) \ (J2*U1);

%%
M = zeros(2*N,n_xhat+1);
b = zeros(2*N,1);
for i = 1:N
    c = Chat*inv(exp(1i*w(i))*eye(n_xhat) - Ahat);
    M(2*i-1,1:n_xhat) = real(c);
    M(2*i  ,1:n_xhat) = imag(c);
    M(2*i-1,n_xhat+1) = 1;
    b(2*i-1)          = real(Gest(i));
    b(2*i  )          = imag(Gest(i));
end

x = M \ b;
Bhat = x(1:n_xhat);
Dhat = x(n_xhat+1:end);

Ghat = tf(ss(Ahat, Bhat, Chat, Dhat, Ts));
% SYS = ss(G);

figure();
semilogy(w, abs(G0)); hold("on");
semilogy(w, abs(squeeze(freqresp(Ghat, w))));
grid();
xlabel("freq [rad/s]")
