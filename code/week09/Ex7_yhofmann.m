%%
%1)
clear all
close all

sigma_u = 1;
sigma_e = 0.1;
number_u = 200;
A = [1 0.25 -0.2 0.1 0.05];
B = [0.6 0.3 -0.05 0];
G = tf(B,A,1);
H = tf(1,A,1);

u_id = sigma_u * randn(number_u,1);
e_id = sigma_e *randn(number_u,1);
v_id = lsim(H,e_id);
y_id = lsim(G,u_id) + v_id;

%%
% 2.)
n_a = 1:8;
n_b = n_a;
system = cell(2,8);


% implement ARX
for i = 1:length(n_a)

    [system{1,i}, system{2,i}] = systemARX(i,i,y_id,u_id);

end


%%
% 3.)
u_val = sigma_u * randn(number_u,1);
e_val = sigma_e * randn(number_u,1);

y_val= lsim(G,u_val) + lsim(H,e_val);

%%
% 4.)

% Predict output for each and plot
y_pred_val = cell(1,8);
for i =1:length(n_a)
    func = tf(system{2,i}',[1,system{1,i}'],1);
    y_pred_val{i} = lsim(func,u_val);
    figure(i)
    plot(y_pred_val{i});
    hold on
    plot(y_val);
    title(sprintf('Predicted for n_a=n_b = %d',i))
    legend('Predicted', 'True');
end


% calculate for identification
y_pred_id = cell(1,8);
for i =1:length(n_a)
    func = tf(system{2,i}',[1,system{1,i}'],1);
    y_pred_id{i} = lsim(func,u_id);
end

%%
% 5.)

% Compute error for each and plot
error_val = zeros(1,8);
error_id = zeros(1,8);
for i = 1:length(n_a)
    error_val(1,i) = mean((y_val - y_pred_val{1,i}).^2);
    error_id(1,i) = mean((y_id - y_pred_id{1,i}).^2); 

end

figure(9)
plot(n_a,error_val);
hold on
plot(n_a, error_id);
xlabel('Order ARX')
title('Error')
legend('Validation', 'Identification')


%%
%6.)

% Compute residuals for each and plot
residual = zeros(number_u,8);
for i = 1: length(n_a)
      residual(:,i) = (y_val - y_pred_val{1,i});
end
figure(100)
plot(1:number_u,residual);
title('Residual for different order')
legend('ARX1', 'ARX2', 'ARX3', 'ARX4', 'ARX5', 'ARX6', 'ARX7', 'ARX8')


%%
%7.)

% Compare freq response
figure (200)
bode(G)
hold on

for i = 1:length(n_b)
    bode(tf(system{2,i}',[1, system{1,i}'],1),u_val)
end
hold off
legend('True', 'ARX1', 'ARX2', 'ARX3', 'ARX4', 'ARX5', 'ARX6', 'ARX7', 'ARX8')

%%
%Functions
function [A,B] = systemARX(orderA, orderB, y, u) % can actually assume that was 0 before so no truncation for u and y needed
phi = zeros(length(y),orderA+orderB);
for i = 1:orderA
    y_temp = -[zeros(i,1);y(1:end-i)];
    phi(:,i) = y_temp;
end

for j = 1:orderB
    u_temp =[zeros(j,1);u(1:end-j)];
    phi(:,orderA + j) = u_temp;
end

theta = phi\y;
A = theta(1:orderA);
B = theta(orderA +1:end);



%could also calculte G here and then make this the return value
%G = tf(B', [1, A'], 1);
end