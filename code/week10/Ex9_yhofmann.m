clear all, close all

N = 10^3;
A = [1 -1.5 0.7];
B = [0 1 0.4];
C = [1 -1.1 0.4];


e = randn(N,1);
u = getU(N,e);
y_true = lsim(tf(B,A,1),u) + lsim(tf(C,A,1),e);  
y_true_filtered = lsim(tf(1,C,1),y_true);
u_filtered = lsim(tf(1,C,1),u);

%%
%1.)
[A_LS, B_LS] = LS_from_ARMAX(N,y_true_filtered,u_filtered);

phi = zeros(N,4);
phi(2,:) = [-y_true_filtered(1) 0 u_filtered(1) 0];
for k = 3:N
    phi(k,:) = [-y_true_filtered(k-1) -y_true_filtered(k-2) u_filtered(k-1) u_filtered(k-2)];
end

theta = phi\y_true_filtered;


%%
%2.)
e_val = randn(N,1);
u_val = getU(N,e_val);
%u_val_filtered = lsim(1/C,u_val);

y_val = lsim(tf(B,A,1),u_val) + lsim(tf(C,A,1),e_val);
y_pred = lsim(tf(B_LS,A_LS,1), u_val) + lsim(tf(C,A_LS,1),e_val);

figure(1)
plot(y_val)
hold on
plot(y_pred)
legend('True', 'Predicted')
xlabel('step k')
ylabel('value')



%%
%3.)

num_Exp = 1000;
A_3 = zeros(3,num_Exp);
B_3 = zeros(3,num_Exp);


for k = 1:num_Exp
    e_exp = randn(N,1);
    u_exp = getU(N,e_exp);
    
    y_exp = lsim(tf(B,A,1),u_exp) + lsim(tf(C,A,1),e_exp);  
    y_exp_filtered = lsim(tf(1,C,1),y_exp);
    u_exp_filtered = lsim(tf(1,C,1),u_exp);

    [A_3(:,k), B_3(:,k)] = LS_from_ARMAX(N,y_exp_filtered,u_exp_filtered);


end


figure(2)
sgtitle('Histograms LS estimates different noise realizations (n = 1000)')

subplot(2,2,1)
histogram(A_3(2,:))
title('a_1')

subplot(2,2,2)
histogram(A_3(3,:))
title('a_2')

subplot(2,2,3)
histogram(B_3(2,:))
title('b_1')


subplot(2,2,4)
histogram(B_3(3,:))
title('b_2')






%%
%FUNCTIONS

function u = getU(N, e)
    u = zeros(N,1);
    u(2) = 0.9*e(1);
    for k = 3:N
        u(k) = 0.14 *u(k-1) +0.12*u(k-2) + 0.9 * e(k-1) + 0.3 * e(k-2);
    end
end


function [A,B] = LS_from_ARMAX(N,y,u)

phi = zeros(N,4);
phi(2,:) = [-y(1) 0 u(1) 0];
for k = 3:N
    phi(k,:) = [-y(k-1) -y(k-2) u(k-1) u(k-2)];
end

theta = phi\y;

B = [0, theta(3) , theta(4)];
A = [1 , theta(1) , theta(2)];

end

