clear all, close all
%% Load and initialize
data = importdata("p3_data.mat");
r1 = data.p3_r1;
y1 = data.p3_y1;

z = tf('z');
F = 0.5*z^-2 + z^-1;
F2 = tf([0, 1,0.5],1,-1, 'Variable','z^-1'); % also like this to avoid a fractionn
%%
%1)
u1 = r1 - lsim(F,y1);

[A,B] = systemARX_no_Initial_Conditions(2,1,y1,u1);
a1 = A(1) 
a2 = A(2)
b1 = B

G_est = tf([b1 0],[1 a1 a2],1);
Est_closed_system = feedback(G_est,F);
y_est = lsim(Est_closed_system,r1);

figure(1)
plot(y1);
hold on
plot(y_est)
title('Direct Method')
legend('Data', 'Estimated')
grid on



%%
% 2.)
% MATLAB code would be exactly the same as above but error term (which does
% not matter for calculation is different and is now be unbiased)


%%
%Functions

%!!!!!! Instead of only writng it for the truncated case can also use the
%normal one and then simply for the \ solving disregard the first 3 rows
function [A,B] = systemARX_no_Initial_Conditions(orderA, orderB, y, u) % here cannot assume that initial conditition was 0 
phi = zeros(length(y) - orderA -orderB+1,orderA+orderB);
for i = 1:orderA
    y_temp = -[y(orderA+orderB - i:end-i)];
    phi(:,i) = y_temp;
end

for j = 1:orderB
    u_temp =[u(orderA + orderB -j:end-j)];
    phi(:,orderA + j) = u_temp;
end

theta = phi\y(orderA+orderB:end);
A = theta(1:orderA);
B = theta(orderA +1:end);

end