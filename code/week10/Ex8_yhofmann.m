clear all, close all
%%
%a)
[u, y] = HS2023_SysID_Exercise_08_GenerateData(20914032);

%%
%b)

r = 20;
g_LSQ = timeEstimation(u,r,y);



%%
%c)
alpha = 0.5;
gamma = 100;

P = kernelTC(r,alpha);
phi = calculatePhi(u,r);
g_hat = (phi'*phi+ gamma*TCKernel_inv(alpha,r))\(phi'*y);



%%
%d)
u_train = u(1:0.7*length(u));
u_test = u(0.7*length(u)+1:end);

y_train = y(1:0.7*length(y));
y_test = y(0.7*length(y)+1:end);


%%
%e)
%TODO
r_candidates = 2:1:30;
alpha_candidates = 0:0.05:0.95;
find_best_alpha_r = zeros(length(alpha_candidates),length(r_candidates),2);
for i = 1:length(alpha_candidates)
    find_best_alpha_r(i,:,1) = r_candidates;
end
for j = 1:length(r_candidates)
    find_best_alpha_r(:,j,2) = alpha_candidates';
end


tuned_alpha = 0;
tuned_r = 2;
tuned_g_hat = cell(1,1);
smallest_error = 10e20;

error_ = zeros(length(alpha_candidates),length(r_candidates));
%g_hat_tuned = cell(length(alpha_candidates),length(r_candidates));
%y_hat_tuned = zeros(length(u_train),size(alpha_candidates) )

for i= 1:length(alpha_candidates)
    for j = 1:length(r_candidates)
    alpha = find_best_alpha_r(i,j,2);
    r = find_best_alpha_r(i,j,1);

    P = kernelTC(r,alpha);
    phi = calculatePhi(u_train,r);
   % g_hat_tuned{i,j} = (phi'*phi+ gamma*inv(P))\(phi'*y);
   g_hat_current = (phi'*phi + gamma*TCKernel_inv(alpha,r))\(phi'*y_train);
   y_estimated = pulseResponseOutput(u_test,g_hat_current);
   error_(i,j) = mean((y_estimated-y_test).^2); % debug
    if(mean((y_estimated-y_test).^2) < smallest_error)
        smallest_error = mean((y_estimated-y_test).^2)
        tuned_r = r;
        tuned_alpha = alpha;
        tuned_g_hat{1,1} = g_hat_current;

        if(smallest_error) < 1
        figure(21436)
        plot(y_estimated)
        hold on
        plot(y_test)
        legend('estimated', 'true')
    i = 2
    close all
    end
    end
    end
end

disp("Best alpha = " + num2str(tuned_alpha));
disp("Best r = " + num2str(tuned_r));





%%
%f)
%TODO

figure(100)
stem(g_LSQ)
hold on
stem(g_hat)
hold on
stem(tuned_g_hat{1,1})

title('Different pulse responses')
xlabel('parameter number')
ylabel('Value')
legend('Least square','Regularized proposed', 'Regularized tuned')



figure(200)
surf(error_)
xlabel('order r')
ylabel('alpha')
zlabel('error')
title('Error for different parameter values')

%% FUNCTIONS
function g = timeEstimation(u, tau_max, y)
    c_toep = u;
    r = [u(1), zeros(1,tau_max)];
    phi = toeplitz(c_toep,r);
    
 

    g = phi\y;
    % other approach
    % phy2 = zeros(length(u), tau_max);
    % for i = 1:tau_max+1
    %     phy2(i:end,i) = u(1:end-i+1);
    % end
    % g2 = phy2\y;

end


function P = kernelTC(N,alpha) 
P = zeros(N,N);
for i = 1:N
    for j = 1:N
        P(i,j) = alpha^max(i,j);
    end
end
end


function Phi = calculatePhi(u,N)
    c_toep = u;
    r = [u(1), zeros(1,N-1)];
    Phi = toeplitz(c_toep,r);



end



function y = pulseResponseOutput(u,pulse)
    tau_max = length(pulse) -1;
    c_toep = u;
    r = [u(1), zeros(1,tau_max)];
    phi = toeplitz(c_toep,r);

    y = phi * pulse;

end


function invP = TCKernel_inv(alpha, n)
diag0 = alpha .^ (-(1:n)') / (1-alpha);
diag0(2:n-1) = diag0(2:n-1) * (1+alpha);
diag1 = -alpha .^ (-(1:n-1)') / (1-alpha);
invP = diag(diag0, 0) + diag(diag1, +1) + diag(diag1, -1);
end