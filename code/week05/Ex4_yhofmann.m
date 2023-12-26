close all
clear all
clc
%%
%a)
N = 1024;
e = randn(N,1)*sqrt(0.01); % pay attention is square root for the standard deviation!!!!!!!!
u = 2 * randi([0, 1], 1, N) - 1;


T = 1;
z = tf('z', T);
H = 0.5*(z-0.9)/(z-0.25)
G = (0.1*z)/(z^4-2.2*z^3+2.42*z^2-1.87*z+0.7225)
time_steps = T*(0:N-1);
resu = lsim(G,u,time_steps);
rese = lsim(H,e,time_steps);
y_experimet = resu + rese;
omega = (2*pi/(T*N)*(0:N-1))';
idx = find(omega>=0 & omega< pi); 

% Uncomment to see plots of signals
%{
figure(1)
plot(u)
hold on
plot(e)
hold on
plot(y_experimet)
legend('u', 'e', 'Y')


figure(2)
plot(rese)
hold on
plot(resu)
hold on
plot(y_experimet)
hold on
legend('H(e)', 'G(u)', 'Y')
%}

%%
%b)
U = fft(u);
Y_exp = fft(y_experimet);

% if one of them has many which are 0 is probably due to periodic input -->
% only calculate at these

 % figure(1)
 % loglog(omega(idx), abs(U(idx)));
% title('U')
% figure(2)
% loglog(omega(idx), abs(Y_exp(idx)));
% title('Y_exp')

G_est = Y_exp ./transpose(U); %Watch out as U is complex don't do U' but use transpose(U)!!!!!!!!!
G_fresp_true = squeeze(freqresp(G,omega));


figure(3)
loglog(omega(idx),abs(G_est(idx)))
hold on
loglog(omega(idx),abs(G_fresp_true(idx)))
hold on
Err = G_est - G_fresp_true;
loglog(omega(idx),abs(Err(idx)))
legend('Estimate', 'True', 'Error');
title('Magnitude')




%%
%c)
N_c = N/4;
omega_c = (2*pi/(T*N_c)*(0:N_c-1))';
idx_c = find(omega_c>=0 & omega_c<= pi);

% Does the same but not really good to do it like that
Y_experimet1 = fft(y_experimet(1:end/4)); 
Y_experimet2 = fft(y_experimet(end/4+1:end/2));
Y_experimet3 = fft(y_experimet(end/2+1:3*end/4));
Y_experimet4 = fft(y_experimet(3*end/4+1:end));
U1 = fft(u(1:end/4));
U2 = fft(u(end/4+1:end/2));
U3 = fft(u(end/2+1:3*end/4));
U4 = fft(u(3*end/4+1:end));
G1_est = Y_experimet1./transpose(U1);
G2_est = Y_experimet2./transpose(U2);
G3_est =Y_experimet3./transpose(U3);
G4_est =Y_experimet4./transpose(U4);

%Better like that
u_res = reshape(u,N_c,4);
y_res = reshape(y_experimet,N_c,4);

U_res = fft(u_res,[],1);
Y_res = fft(y_res,[],1);

G_est_res = Y_res./U_res;


G_est_average = zeros(length(G1_est),1);
for k =1:length(G1_est)
    G_est_average(k) = mean([G1_est(k);G2_est(k);G3_est(k); G4_est(k)]);
end
G_sum = G1_est+G2_est+G3_est+G4_est;
G_est_average2 = (G_sum)/4; % both ways work now


G_est_average2 = mean(G_est_res,2);


figure(5)
loglog(omega_c(idx_c),abs(G_est_average2(idx_c)))
hold on
loglog(omega(idx),abs(G_est(idx)))
hold on
loglog(omega(idx),abs(G_fresp_true(idx)))
legend('Averaged estimate', 'Estimate', 'True System');
title('Magnitude from averaged')
